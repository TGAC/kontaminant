

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
#include <locale.h>
#include <assert.h>
#include <limits.h>
#ifdef THREADS
#include <pthread.h>
#endif

#include <global.h>
#include <flags.h>
#include <nucleotide.h>
#include <binary_kmer.h>
#include <seq.h>
#include <element.h>
#include <open_hash/hash_table.h>
#include <dB_graph.h>
#include <path.h>
#include <perfect_path.h>
#include <y_walk.h>
#include <branches.h>
//#include <binary_tree.h>
#include <mark_pair.H>
#include <logger.h>
static ReadConnectionGraph * mark_pair_get_links(int library,  dBGraph * db_graph);

int mark_pair_asign_perfect_paths(Path * path, dBGraph * db_graph, uint32_t id){
  
    pathStep current;
    int marked = 0;
    int i;
    path_get_last_step(&current, path);
    if (db_node_edges_count(current.node, forward) + db_node_edges_count(current.node, forward) == 1) {
        db_node_set_supernode_id(id, current.node);
        marked++;
    }
    path_get_step_at_index(0, &current, path);
    if (db_node_edges_count(current.node, forward) + db_node_edges_count(current.node, forward) == 1) {
        db_node_set_supernode_id(id, current.node);
        marked++;
    }  
    for(i = 1; i < path_get_length(path)-1; i++){
        path_get_step_at_index(i, &current, path);
        db_node_set_supernode_id(id, current.node);
        marked++;
    }
    
    return marked; 
}

static void mark_adjacency_from_node(dBNode * node, void * argsa){
    mark_args * args = (mark_args *)argsa;
    Path * path = NULL;
    int len;
    int mark;
    if (db_node_get_supernode(node) == 0) {
        path = path_get_buffer_path();
        len = perfect_path_get_path(node, undefined, &db_node_action_do_nothing, args->graph, path );
        if(len > 0){
            mark = mark_pair_asign_perfect_paths(path, args->graph, args->count +1);
            if(mark)
                args->count++;    
            else
                args->fail_count++;
        }
        path_free_buffer_path(path);
        path = NULL;
    }    
   
}

long long mark_pair_mark_contiguous_perfect_paths(dBGraph * db_graph){
    int threads =  1, i;
#ifdef THREADS
    threads = db_graph->number_of_threads;
#endif
    long long total_count = 0;
    long long total_fail = 0;
    mark_args ** args = calloc(threads, sizeof(mark_args  *));
    for (i=0; i<threads; i++) {
        args[i] = calloc(1, sizeof(mark_args ));
        args[i]->graph = db_graph;
        args[i]->count = 0;
        args[i]->fail_count = 0;
    }
    log_write_timestamp(true);
    log_and_screen_printf("Marking perfect paths\n");
    hash_table_traverse_with_args(&mark_adjacency_from_node, (void **)args, db_graph);
    
    for (i = 0; i < threads; i++) {
        total_count +=  args[i]->count;
        total_fail  +=  args[i]->fail_count;
        free(args[i]);
    }
    
    db_graph->perfect_path_count = total_count;
    log_and_screen_printf("Marked %lld perfect paths, not marked %lld perfect paths\n", total_count, total_fail);
    free(args);
    return total_count;

}

void mark_pair_validate_file(char * filename){
    FILE *fp_fof = fopen(filename, "r");
    FILE *fp_1;
    FILE *fp_2;
    char file_1[1000];
    char file_2[1000];
    int distance, library;
    
    if (!fp_fof) {
        printf("Error: can't open file %s\n", filename);
        exit(-1);
    }
    
    while (!feof(fp_fof)) {
        fscanf(fp_fof, "%s %s %d %d\n", file_1, file_2, &distance, &library);
        
        if (strcmp(file_1, file_2) == 0) {
            log_and_screen_printf("Error: Both files in the pair are the same!\n");
            fclose(fp_fof);
            exit(-1);
        }
        
        fp_1 = fopen(file_1, "r");
        if (!fp_1) {
            log_and_screen_printf("Error: can't open file %s\n", file_1);
            fclose(fp_fof);
            exit(-1);
        }
        
        fp_2 = fopen(file_2, "r");
        if (!fp_2) {
            log_and_screen_printf("Error: can't open file %s\n", file_2);
            fclose(fp_fof);
            fclose(fp_1);
            exit(-1);
        }
        
        fflush(stdout);
        
        fclose(fp_1);
        fclose(fp_2);        
    }

}

void mark_pair_load_pair_data(char* filename, char fastq_ascii_offset,dBGraph * db_graph) {
    FILE *fp_fof = fopen(filename, "r");
    
    char file_1[1000];
    char file_2[1000];
    int library, distance;    
    mark_pair_validate_file(filename);    
    
#ifdef THREADS
    hash_table_threaded_traverse(&db_node_action_clear_flags, db_graph);
#else
    hash_table_traverse(&db_node_action_clear_flags, db_graph);
#endif
    
    mark_pair_mark_contiguous_perfect_paths(db_graph);    
    rewind(fp_fof);
    
    while (!feof(fp_fof)) {
        fscanf(fp_fof, "%s %s %d %d\n", file_1, file_2, &distance, &library );        
        printf("Marking pair: %s and %s (distance: %d, library: %d)\n", file_1, file_2, distance, library);
        
        mark_pair_search_enrich_graph(library, file_1, file_2,fastq_ascii_offset,db_graph);
        ReadConnectionGraph * rcg = mark_pair_get_links(library, db_graph);
        rcg->expected_distance = distance;
               
    }
    
    fclose(fp_fof);
}

static int binary_kmer_sliding_window_to_path(KmerSlidingWindowSet * windows, Path * path,dBGraph * db_graph ){
    path_reset(path);
    int count = 0;
    dBNode * current_node;
	pathStep ps;
	ps.label = Undefined;
	ps.orientation = forward;
    ps.flags = 0;
    BinaryKmer tmp_key;
    short kmer_size = db_graph->kmer_size;
    
    while (binary_kmer_sliding_window_set_get_next(windows)) {
        element_get_key(&windows->current,kmer_size, &tmp_key); 		
        current_node = hash_table_find(&tmp_key, db_graph);
        if (current_node != NULL) {
            ps.node = current_node;
            if (binary_kmer_comparison_operator(tmp_key,windows->current)) {
                ps.orientation = forward;
            } else {
                ps.orientation = reverse;
            }
            path_add_node(&ps, path);
            count++;
        } 
    }
    return count;
}

static void * new_connection_from(void * c){
    PathConnection * connection = c;
    PathConnection * conn = calloc(1,sizeof(PathConnection));
    conn->to = connection->to;
    conn->from = connection->from;
    conn->coverage = connection->coverage;
    return conn;
}

static int compare_connections(void * va, void * vb){
    PathConnection * a = va;
    PathConnection *b = vb; 
    if (a->to == b->to && a->from == b->from) {
        return 0;
    }
    
    if (a->to < b->to) {
        return  -1;
    }
    if (a->to > b->to) {
        return  1;
    }
    
    if (a->from < b->from) {
        return  -1;
    }
    if (a->from > b->from) {
        return  1;
    }
    fprintf(stderr, "The comparasion is invalid.");
    assert(0);
    
}

/*static BinaryTree * mark_pair_get_links(dBGraph * db_graph){
    if (db_graph->supernode_links == NULL) {
        db_graph->supernode_links = binary_tree_new(db_graph->bucket_size, &compare_connections);//The minumum shall be calculated in better way? 
        
    }
    return db_graph->supernode_links;
}*/

static unsigned long db_graph_get_perfect_path_count(dBGraph * db_graph ){
    assert(db_graph->perfect_path_count != 0);
    return db_graph->perfect_path_count;
}

static void db_graph_init_libraries(dBGraph * db_graph){
    if (db_graph->supernode_link == NULL) {
        db_graph->supernode_link = calloc(MAX_READ_LIBRARIES, sizeof(ReadConnectionGraph));
        if (db_graph->supernode_link == NULL) {
            fprintf(stderr, "Unable to alloc supernode_link\n");
            exit(-1);
        }
    }
}

static void init_connection_graph(int library, dBGraph * db_graph){
    unsigned long path_count = db_graph_get_perfect_path_count(db_graph);

    ReadConnectionGraph * rcg = db_graph->supernode_link ;
    
    if( rcg[library].links == NULL) {
        rcg[library].links = calloc(path_count *  path_count, sizeof(mp_cov_counter));
        rcg[library].number_of_paths = path_count;
    }
    
}


//TODO: make a multi library implementation. At the moment, there is only one. 
static ReadConnectionGraph * mark_pair_get_links(int library,  dBGraph * db_graph){
    if (db_graph->supernode_link == NULL) {
        db_graph_init_libraries(db_graph);
    }
    ReadConnectionGraph * rcg = db_graph->supernode_link ;
    if (rcg[library].links == NULL) {
        init_connection_graph(library, db_graph);
    }
    return &rcg[library];
}

static mp_cov_counter mark_pair_connect(int library, uint32_t p1, uint32_t p2, dBGraph * db_graph){
    if (p1 == p2) {
        return  0;//We havent added a new link. 
    }
    
    ReadConnectionGraph * links = mark_pair_get_links(library, db_graph);
    assert(links != NULL);
    
    return mark_pair_connections_between(p1, p2, links);   

}

static int mark_pair_distance_between_pairs(Path * p1, Path * p2, dBGraph * db_graph){
    pathStep first, second;
    path_get_step_at_index(0, &first, p1);
    path_get_last_step(&second, p2);
    
    Path * buff = path_get_buffer_path();
    if (first.node->supernode == second.node->supernode) {
        perfect_path_get_path(first.node, first.orientation, &db_node_action_do_nothing, db_graph, buff);
        return path_index_of_step(&second,buff ); //TODO: What with the opposite orientation? maybe we
    }
    
    path_free_buffer_path(buff); 
    return -1;
    
}

static int  mark_pair_connect_supernodes_between(int library, Path * p1, Path * p2, dBGraph * db_graph){
    assert(p1 != NULL);
    assert(p2 != NULL);
    
    int count = 0;
    int i, j;
    ReadConnectionGraph * rcg = mark_pair_get_links(library, db_graph);
    
    for (i=0; i <p1->supernodes_total_count; i++) {
        for (j=0; j < p2->supernodes_total_count; j++) {	
            if (mark_pair_connect(library, p1->supernodes[i], p2->supernodes[j],db_graph)) {
                count++;
            }
        }
    }
    if (count) {
        rcg->used_reads++;
    }else{
        rcg->unused_reads++;
    }
    return count;
}

/*
uint32_t mark_pair_connections_between(uint32_t p1, uint32_t p2, dBGraph * db_graph){
    assert(p1 != p2);    
    uint32_t count = 0;
    
    PathConnection key;
    key.from = p1;
    key.to = p2;
    
    PathConnection * link = (PathConnection * )binary_tree_find(&key, db_graph->supernode_links);
    if (link != NULL) {
        count = link ->coverage;
    }
    
    return count;
}*/

mp_cov_counter mark_pair_connections_between(uint32_t p1, uint32_t p2, ReadConnectionGraph * rcg){
    assert(p1 != p2);    
    mp_cov_counter count = 0;
    uint32_t shift_size = rcg->number_of_paths;
    uint32_t offset = p1 * shift_size + p2;
    count = rcg->links[offset];
    if (rcg->links[offset] == 0) {
        rcg->connected_paths++;
    }
    if ( count < MAX_MP_COV) {
        count = rcg->links[offset]++;
    }else{
        rcg->overflows++;
    }
    return count;
}

void mark_pair_log_connection_stats(dBGraph * db_graph){
    int i;
    for (i = 0; i < MAX_READ_LIBRARIES; i++) {
        ReadConnectionGraph * rcg = mark_pair_get_links(i, db_graph);
        if (rcg != NULL) {
            log_and_screen_printf("Library %d\n", i);
            log_and_screen_printf("Number of paths %d\n", rcg->number_of_paths);
            if(rcg->distance_sum){
                log_and_screen_printf("Average distance %d\n", rcg->distance_sum/rcg->used_reads);
            }
            log_and_screen_printf("Used reads %d\n", rcg->used_reads);
            log_and_screen_printf("Connected_paths %d (%f%%)\n", rcg->connected_paths, 100*(float)rcg->connected_paths/((float)rcg->number_of_paths*(float)rcg->number_of_paths));
            log_and_screen_printf("Expected distance %d\n", rcg->expected_distance);
            log_and_screen_printf("Overflows %d\n", rcg->overflows);
            log_and_screen_printf("Distance search range (%d, %d)\n\n", rcg->min_distance, rcg->max_distance);

        }        
    }

}

static int binary_kmer_to_path(KmerSlidingWindowSet * windows, Path * path_read, dBGraph * db_graph){
    short kmer_size = db_graph->kmer_size;
    dBNode * current_node;
    pathStep ps;
    ps.label = Undefined;
    BinaryKmer tmp_key;
    int missed_kmers = 0;
    path_reset(path_read);
    while (binary_kmer_sliding_window_set_get_next(windows)) {
        element_get_key(&windows->current,kmer_size, &tmp_key); 		
        current_node = hash_table_find(&tmp_key, db_graph);
        if (current_node != NULL) {
            ps.node = current_node;
            if (binary_kmer_comparison_operator(tmp_key,windows->current)) {
                ps.orientation = forward;
            } else {
                ps.orientation = reverse;
            }
            path_add_node(&ps, path_read);
        } else {
            missed_kmers++;
        }
    }
    return missed_kmers;
}
//  mark_pair_search_enrich_graph(fp_1, fp_2,fastq_ascii_offset,db_graph);
void mark_pair_search_enrich_graph(int library, char * file_1, char * file_2,  char fastq_ascii_offset, dBGraph * db_graph)
{
    
	FILE * f1;
    FILE * f2;
	int length_1, length_2, nkmers1, nkmers2;
    int missed_kmers_1, missed_kmers_2;
	short kmer_size = db_graph->kmer_size;
    
    
	int max_windows = MAX_READ_LENGTH / (kmer_size + 1);
	
	//number of possible kmers in a 'perfect' read
	int max_kmers = MAX_READ_LENGTH - kmer_size + 1;
	
	KmerSlidingWindowSet * windows = binary_kmer_sliding_window_set_new(max_windows, max_kmers);
	windows->kmer_size = kmer_size;
	
	int quality_cut_off = 0; //The cutoff comes from the previous CTXs, or the cleaning algorithm. 
	Path * path_read_one = path_get_buffer_path();
	Path * path_read_two = path_get_buffer_path();
	
	pathStep ps;
	ps.label = Undefined;
	ps.orientation = forward;
    ps.flags = 0;
    long long seq_count = 0;
	
    long long marked_on_same_read = 0;
    long long marked_on_different_read = 0;
    
    Sequence * seq1 = sequence_new(MAX_READ_LENGTH, MAX_FILENAME_LENGTH, fastq_ascii_offset);
	Sequence * seq2 = sequence_new(MAX_READ_LENGTH, MAX_FILENAME_LENGTH, fastq_ascii_offset);
    
    f1 = fopen(file_1, "r");
    if (!f1) {
        log_and_screen_printf("Error: can't open file %s\n", file_1);
        exit(-1);
    }
    
    f2 = fopen(file_2, "r");
    if (!f2) {
        log_and_screen_printf("Error: can't open file %s\n", file_2);
        fclose(f1);
        exit(-1);
    }

    int max_read_length = 0;
	while((length_1 = read_sequence_from_fastq(f1, seq1, MAX_READ_LENGTH) )&& (length_2 = read_sequence_from_fastq(f2, seq2, MAX_READ_LENGTH))){
		seq_count++;
        path_reset(path_read_one);
        path_reset(path_read_two);
        
        if ((seq_count % 1000000) == 0) {
            fprintf(stdout, ".");
        }
        
		
        missed_kmers_1 = 0;
        nkmers1 = get_sliding_windows_from_sequence(seq1->seq, seq1->qual, seq1->length, quality_cut_off, kmer_size, windows, max_windows, max_kmers,false,0);					
		binary_kmer_to_path(windows, path_read_one, db_graph);
        
        missed_kmers_2 = 0;
        nkmers2 = get_sliding_windows_from_sequence(seq2->seq, seq2->qual, seq2->length, quality_cut_off, kmer_size, windows, max_windows, max_kmers, false, 0);
        binary_kmer_to_path(windows, path_read_two, db_graph);
        
        //We connect between reads and on the same read. 
        marked_on_different_read += mark_pair_connect_supernodes_between(library,path_read_one,path_read_two, db_graph);
        if (library == 1) {
            if (max_read_length<nkmers1) {
                max_read_length = nkmers1;
            }
            if (max_read_length < nkmers2) {
                max_read_length = nkmers2;
            }
            marked_on_same_read += mark_pair_connect_supernodes_between(0,path_read_one, path_read_one, db_graph);
            
            marked_on_same_read += mark_pair_connect_supernodes_between(0, path_read_two, path_read_two, db_graph);
        }
       
	}
    
    ReadConnectionGraph * rcg_0 = mark_pair_get_links(0, db_graph);
    if (rcg_0->expected_distance < max_read_length) {
        rcg_0->expected_distance = max_read_length;
    }
    
   	path_free_buffer_path(path_read_one);
    path_free_buffer_path(path_read_two);
   
    
	binary_kmer_free_kmers_set(&windows);
	free_sequence(&seq1);
	free_sequence(&seq2);
    
    
    log_write_timestamp(1);
    log_and_screen_printf("\nEnriched %lli sequences\n ", seq_count);
    log_and_screen_printf("\t Marked:\t%lli \n ", marked_on_different_read + marked_on_same_read);
    log_and_screen_printf("\t same read:\t %lli \n ", marked_on_same_read);
    log_and_screen_printf("\t different read:\t%lli \n ", marked_on_different_read);
    
  
}


/**
*
* Here starts the code for the actual walk. 
*
*
*/

static void pre_step_action(pathStep * ps){
	
	if(ps->orientation == forward){
        db_node_action_set_flag(ps->node, VISITED_FORWARD);
	}else{
        db_node_action_set_flag(ps->node, VISITED_REVERSE);
	}
	
}

static void post_step_action(pathStep * ps){
	if(ps->orientation == forward){
        db_node_action_unset_flag(ps->node, VISITED_FORWARD);
		
	}else{
        db_node_action_unset_flag(ps->node, VISITED_REVERSE);
	}
}

//TODO: encapsulate this...
static long long count_solved_x_nodes = 0;
static long long count_unsolved_x_nodes = 0;
static long long count_solved_short_x_path = 0;
static long long count_unsolved_short_x_path = 0;
static long long count_solved_long_x_path = 0;
static long long count_unsolved_long_x_path = 0;

static void increas_global_counter(long long * counter){//Wrapper to make this multithreaded. 

    *counter += 1;
}


static Nucleotide mark_pair_get_best_next_step_nucleotide_for_X_node(dBNode * from, dBNode * previous,  Orientation orientation, dBGraph * db_graph ){
    assert(from != NULL);
    assert(previous != NULL);
    Nucleotide best = Undefined;
    
    dBNode * neighbours[4];
    Nucleotide labels[4];
    int found = db_graph_get_neighbouring_nodes_all_colours(from, orientation, neighbours, labels, db_graph);
    int i;
    
    uint32_t conn_count = 0 ;
    uint32_t curr_count;
    boolean valid = true;
    uint32_t supernode_prev = db_node_get_supernode(previous);
    uint32_t supernode_curr = 0;
    if (supernode_prev == 0) {
        valid = false;
    }
    for (i = 0; i < found && valid; i++) {
        supernode_curr = db_node_get_supernode(neighbours[i]);                                            
        if (supernode_curr != 0) {
            curr_count = mark_pair_connections_between(supernode_prev, supernode_curr, mark_pair_get_links(0, db_graph));//Note: Library 0 contains the nodes connected from reads. 
            if (  conn_count && curr_count) {
                valid = false;  
            }
            if (conn_count == 0 && curr_count ) {
                conn_count = curr_count;
                best = labels[i];
            }
        }
    }
    if (valid == false) {
        best = Undefined;
    }
       
    return best;
}

static Nucleotide mark_pair_get_best_nucleotide_from_history(int library, pathStep * current_step, pathStep * next_step, dBGraph * db_graph){
    Path * hist = current_step->path;
    Nucleotide best = Undefined;
    pathStep hist_step;
    
    ReadConnectionGraph * rcg = mark_pair_get_links(library, db_graph);
    int max_distance = rcg->expected_distance;
    
    int hist_length = path_get_length(hist);
    int wall = hist_length - max_distance;
    wall = wall < 0? 0:wall;
    int last_index = path_get_first_in_node_after(wall, hist);
    int i;
       if (last_index != -1) {
        int  best_matches = 0, temp_matches, second_best_matches = 0;
        int historic_path, current_path = current_step->node->supernode;
        
        if (last_index > 0) {
            do {
                last_index--;
                path_get_step_at_index(last_index, &hist_step, hist);
                historic_path = hist_step.node->supernode;
            } while (historic_path == 0 && last_index > wall);
            
        }
        
        if (last_index > wall) {
            dBNode * neighbours[4];
            Nucleotide labels[4];
            int found = db_graph_get_neighbouring_nodes_all_colours(next_step->node, next_step->orientation, neighbours, labels, db_graph);
            boolean invalid; //Maybe we want another counter for this case, when we see that the the double Y is connecting to the next path. 
            
            for (i = 0; i < found; i++) {
                current_path = neighbours[i]->supernode;
                if (current_path == historic_path) {
                    invalid = true;
                }else{
                    temp_matches = mark_pair_connections_between(current_path, historic_path, rcg);
                    if (temp_matches) {
                        second_best_matches = best_matches;
                        best_matches = temp_matches;
                        best = labels[i];
                    }
                }
            }
            
            
            if (best_matches <= (second_best_matches * 3) || invalid) {//TODO: This shall be a variable, of how much difference we want between both paths.
                best_matches = Undefined;
            }
        }
        

    }
        
    return best;
}

static pathStep *get_next_step(pathStep * current_step, pathStep * next_step,
                               pathStep * reverse_step, dBGraph * db_graph)
{
	pathStep *step = db_graph_get_next_step(current_step, next_step, reverse_step,
                                            db_graph);
	assert(step != NULL);
    Nucleotide n;
	int count_fwd, count_rev;
    if (step->node != NULL) {
        
		step->label = Undefined;
        count_fwd = db_node_edges_count(next_step->node, next_step->orientation);
        count_rev = db_node_edges_count(next_step->node, opposite_orientation(next_step->orientation));
		if(count_fwd == 0){
            step->label = Undefined;//We cant do anything. 
        }else if (db_node_has_precisely_one_edge(next_step->node, next_step->orientation, &n)) {//Simple case, no ambiguity, we don't need to do anything else. 
			next_step->label = n;
		}else if(count_fwd > 1 && count_rev > 1 ){ //This is an X node, try to resolve with the own reads tusignature. If it is not greater than 1, we need to get the previous
            //   path_get_step_at_index(path_get_length(current_path) - 2, &prevStep, current_path);
            next_step->label = mark_pair_get_best_next_step_nucleotide_for_X_node(next_step->node, current_step->node, next_step->orientation, db_graph);
            if(next_step->label != Undefined){
                count_solved_x_nodes++;
            }else{
                count_unsolved_x_nodes++;
            }
        }
        //TODO: Make this with a dynamic/user defined distance and for different distances
        if( step->label == Undefined && count_fwd > 1){
            next_step->label = mark_pair_get_best_nucleotide_from_history(0,    current_step,   next_step, db_graph);
            if(next_step->label != Undefined){
                count_solved_short_x_path++;
            }else{
                count_unsolved_short_x_path++;
            }
        }
        if( step->label == Undefined && count_fwd > 1){
            next_step->label = mark_pair_get_best_nucleotide_from_history(1,  current_step,   next_step, db_graph);
            if(next_step->label != Undefined){
                count_solved_long_x_path++;
            }else{
                count_unsolved_long_x_path++;
            }
        }
        
        /*
        if( step->label == Undefined && count_fwd > 1){//Here we try to do the read pair magic. 
            next_step->label = find_best_nucleatide_in_double_y(current_step, next_step, reverse_step,db_graph);
            if(next_step->label != Undefined){
                count_solved_short_x_path++;
            }else{
                count_unsolved_short_x_path++;
            }
        }//This are all the cases. */
        
	}
    
	return next_step;
}

static boolean continue_traversing(pathStep * current_step,
                                   pathStep * next_step,
                                   pathStep * reverse_step, Path * temp_path,
                                   dBGraph * db_graph)
{
    //pathStep first;
    
    boolean cont = true;
    
    boolean next_visited = 
    (db_node_check_for_flag(next_step->node, VISITED_FORWARD) == true && next_step->orientation==forward )||
    (db_node_check_for_flag(next_step->node, VISITED_REVERSE) == true && next_step->orientation==reverse); //TODO: make this walk a bit forwad to see if the read pair is able to solve the little loop
    
    if (next_visited && next_step->label != Undefined) {
        path_action_set_flag(temp_path, IS_CYCLE);
        path_add_stop_reason(LAST, PATH_FLAG_IS_CYCLE, temp_path);
        cont = false;
    }
    
    if (current_step->label == Undefined){
        path_add_stop_reason(LAST, PATH_FLAG_STOP_BLUNT_END, temp_path);
        cont = false;
    }
    
    if(!path_has_space(temp_path)){
        path_add_stop_reason(LAST, PATH_FLAG_LONGER_THAN_BUFFER, temp_path);
        cont = false;
    }
    return cont;
}


static void path_callback_copy(Path * path, void * arg){
    path_copy(arg, path);
}

static void path_callback_append(Path * path, void * arg){
    path_append(arg, path);
}

Path * mark_pair_get_path_single_walk(dBNode * node, Orientation orientation,
                                      void (*node_action) (dBNode * node),
                                      dBGraph * db_graph, boolean both_directions,Path * path)
{
    
    Path *buff1 = path_get_buffer_path();
    assert(buff1 != path);
    path_reset(buff1);
    path_reset(path);	//Make sure that the paths are clean
    
    pathStep first;
    boolean only_one_edge;
    first.node = node;
    first.orientation =  orientation;
    
    WalkingFunctions wf;
    perfect_path_get_funtions(&wf);
    assert(buff1 != path);

    
    db_graph_add_node_action(&wf, node_action);
    wf.continue_traversing = &continue_traversing;
    
    wf.get_next_step    =  &get_next_step;
    wf.pre_step_action  = &pre_step_action;
    wf.post_step_action = &post_step_action;
    
    
    if(orientation != undefined){
        only_one_edge = db_node_has_precisely_one_edge(first.node, first.orientation, &first.label);
        wf.get_starting_step = &get_first_step_identity;
        
    }else {
        first.orientation = reverse;
        only_one_edge = db_node_has_precisely_one_edge(first.node, first.orientation, &first.label);
    }
    
    //we got the first side. We use path, since from here we will append
    //the second half of the path. 
    
    if (both_directions == true) {
        db_graph_add_path_callback_with_args(&wf, (void (*)(Path *, void *)) &path_reverse, path);
    }else{
       db_graph_add_path_callback_with_args(&wf,  &path_callback_copy, path); 
    }
    
    assert(buff1 != path);
    path_reset(buff1);
    
    if (db_node_edges_count(first.node, first.orientation) > 0) {	
        assert(buff1 != path);
        db_graph_generic_walk(&first, buff1, &wf, db_graph);
        assert(buff1 != path);
    } else {
        if (both_directions == true) {
            first.label = Undefined;
            first.orientation = opposite_orientation(first.orientation);
            path_add_node(&first, path); //In this way, we can use the last node samelesly for the second round, without checking if it was able to walk in the first place
        }
    }
    
    if (both_directions == true){		
        //The second time, we just need to append
        
        boolean removed = db_graph_remove_path_callback(&wf, path_reverse );
        assert(removed != 0);
        
        db_graph_add_path_callback_with_args(&wf, &path_callback_append, path);
        
        path_get_last_step(&first, path);//Starting over from the last step in the path. 
        db_node_has_precisely_one_edge(first.node, first.orientation, &first.label);
        
        if(db_node_edges_count(first.node, first.orientation) > 0 ){
            db_graph_generic_walk(&first, buff1, &wf, db_graph);
        }
        
    }
    
    path_free_buffer_path(buff1);
    
    return path;	
}



typedef struct{
    long long count_kmers;
    dBGraph * db_graph;
    Path * path;
    FILE * fout_cov;
    FILE * fout_viz;
    FILE * fout;
    int count_nodes;
    int max_length;
    int count_repetitive;
    int count_sing;
} makr_pair_print_path_args;

static void mark_pair_print_path(dBNode * node, void * args)
{
    int min_fresh_nodes_to_print = 90;//This goes in percentage. 
    makr_pair_print_path_args * mpppa = (makr_pair_print_path_args *) args;
    mpppa->count_kmers++;
    if (db_node_check_flag_visited(node) == false) {
        db_node_action_set_flag_visited(node);
        Path *  found_path = mark_pair_get_path_single_walk(node, undefined, &db_node_action_do_nothing, mpppa->db_graph, true, mpppa->path);
        if (found_path  ) {
            //pathStep last_step;
            //boolean found_next_path = false;
            path_mark_as_visited(found_path);
            if (path_get_nodes_count(mpppa->path) > 1 && path_percentage_new_nodes(found_path) > min_fresh_nodes_to_print) {
                if (mpppa->fout_cov != NULL) {
                    path_to_coverage(mpppa->path, mpppa->fout_cov);
                }
                if (mpppa->fout_viz != NULL) {
                    path_graphviz_line(mpppa->fout_viz, mpppa->path);
                }
                path_to_fasta(mpppa->path, mpppa->fout);
                if (mpppa->path->length == mpppa->max_length) {
                    printf("contig length equals max length [%i] for node_%i\n",
                           mpppa->max_length, mpppa->count_nodes);
                }
                mpppa->count_nodes++;
                path_increase_id(mpppa->path);
            } if(path_percentage_new_nodes(found_path) <= min_fresh_nodes_to_print){
                mpppa->count_repetitive++;
            }else {
                mpppa->count_sing++;
            }
        }
    }
}

void mark_pair_print_paths(char *filename, int max_length, boolean with_coverages, boolean with_viz, dBGraph * db_graph)
{	
    mark_pair_log_connection_stats(db_graph);
    
    makr_pair_print_path_args mpppa;
	Path *path = path_get_buffer_path();	//We will try to use only this buffer path.
	Path *next_path = path_get_buffer_path();
    
	path_reset(path);
	
    
    mpppa.path = path;
    mpppa.count_kmers = 0;
    mpppa.count_nodes = 0;
    mpppa.count_repetitive = 0;
    mpppa.count_sing = 0;
    mpppa.fout = NULL;    
	mpppa.fout_cov = NULL;
	mpppa.fout_viz = NULL;
	mpppa.db_graph = db_graph;
    
	mpppa.fout = fopen(filename, "w");
	
	if (with_coverages) {
		char filename_cov[strlen(filename) + 10];
		sprintf(filename_cov, "%s_cov", filename);
		mpppa.fout_cov = fopen(filename_cov, "w");
	}
	
	if (with_viz) {
		char filename_cov[strlen(filename) + 10];
		sprintf(filename_cov, "%s.viz", filename);
		mpppa.fout_viz = fopen(filename_cov, "w");
		path_graphviz_open_header(mpppa.fout_viz);
	}
	
	mark_double_y(db_graph);	
	log_and_screen_printf("Printing paths\n");
    void * args[1];
    args[0] = &mpppa;
    
	hash_table_traverse_with_args(&mark_pair_print_path, args, db_graph);
	log_and_screen_printf("%'d nodes visited [%'qd singletons, %'qd repetitive]\n", mpppa.count_nodes, mpppa.count_sing, mpppa.count_repetitive);
    log_and_screen_printf("solved x-nodes %'lld\n unsolved x-nodes: %'lld\n", count_solved_x_nodes, count_unsolved_x_nodes);
    log_and_screen_printf("Short supernodes  %'lld\n unsolved nodes: %'lld\n", count_solved_short_x_path, count_unsolved_short_x_path);
    log_and_screen_printf("Long  supernodes %'lld\n unsolved nodes: %'lld\n", count_solved_long_x_path, count_unsolved_long_x_path);
    
	//path_destroy(path);
	fclose(mpppa.fout);
	if (with_coverages) {
		fclose(mpppa.fout_cov);
	}
	if (with_viz) {
		path_graphviz_close_header(mpppa.fout_viz);
		fclose(mpppa.fout_viz);
	}
    
    path_free_buffer_path(path);
    //path_free_buffer_path(buff);
    path_free_buffer_path(next_path);
    
	//	y_node_inited = false;	//Since we can clean the flags before running the node, we need to make sure next time we will find them	
}

