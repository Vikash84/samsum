#include "helper.h"
using namespace std;


vector<std::string> format_matches_for_service(vector<MATCH*> &all_reads, char* &index) {
    /* Parameters:
      * all_reads: A vector populated by MATCH instances that contain information to be converted into strings
     * Functionality:
      * Iterates through all MATCH instances in all_reads and generates two strings for each: the name of the query
      and the relevant fields that are to be returned into Python.
    */
    vector<std::string> query_info;
    char buf[1000];

    for ( vector<MATCH*>::iterator it = all_reads.begin(); it != all_reads.end(); it++)  {
        if ( (*it)->parity > 1 )
            cout << (*it)->query << " " <<  (*it)->start << " " << (*it)->cigar.c_str() << " " << (*it)->parity << " " << (*it)->w << endl;
        if ( strcmp(index, "q") == 0 ) {
            query_info.push_back((*it)->query);
            sprintf(buf, "%s\t%d\t%s\t%d\t%f", (*it)->subject.c_str(), (*it)->start, (*it)->cigar.c_str(), (*it)->parity, (*it)->w);
        }
        if ( strcmp(index, "r") == 0 ) {
            query_info.push_back((*it)->subject);
            sprintf(buf, "%s\t%d\t%s\t%d\t%f", (*it)->query.c_str(), (*it)->start, (*it)->cigar.c_str(), (*it)->parity, (*it)->w);
        }
        query_info.push_back(buf);
    }
     while ( !all_reads.empty()) {
        delete all_reads.back();
        all_reads.pop_back();
     }

    return query_info;
}


void remove_low_quality_matches(vector<MATCH*> &mapped_reads, unsigned int min_map_qual, float &unmapped_weight_sum) {
    /* Parameters:
      * mapped_reads: A vector of MATCH instances that is to be filtered
      * min_map_qual: An integer representing the minimum mapping quality for an alignment to be included
      * unmapped_weight_sum: Reference to a float that tracks the sum weight of fragments
     * Functionality:
      * A new filtered_matches vector is created and all MATCH instances that pass the min_map_qual are appended to it.
      * Their weights are added to the unmapped_weight_sum float since these are no longer returned as MATCHes
    */
    int pre = mapped_reads.size();
    vector<MATCH*> filtered_matches;
    vector<MATCH*>::iterator it = mapped_reads.begin();
    while (it != mapped_reads.end()) {
        if ((*it)->mq < min_map_qual) {
            unmapped_weight_sum += (*it)->w;
            it = mapped_reads.erase(it);
            filtered_matches.push_back(*it);
        }
        else
            it++;
    }
    mapped_reads.shrink_to_fit();

//    cout << pre << " " << mapped_reads.size() << endl;
//    while ( !filtered_matches.empty()) {
//        delete filtered_matches.back();
//        filtered_matches.pop_back();
//    }

    return;
}
