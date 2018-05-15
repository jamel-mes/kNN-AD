#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <vector>
#include <boost/tokenizer.hpp>
#include <boost/regex.hpp> 
#include "AD_kNN.hpp"

using namespace std ;
using namespace boost ;

#define MOLIDPOS 0
#define DESCRIPTORPOS 1
#define SHRINK 1.3

//Formula : D=<Di>+Z*sigma

int main(int argc, char *argv[])
{
    // Verify if arguments are OK
    if(argc != 5){
        sethelp(argv[0]) ;
        exit(EXIT_FAILURE);
    }

    int K = atoi(argv[3]) ; // K nearest neighbors
	double Z = atof(argv[4]) ; // Z empirical param
/*
 1   Read Training Set file
 */
    vector<string> TR_molid                ; // training set mol identifier
    vector<vector<double> > TR_descriptors ; // training set descriptors
    string line ; int n = 0 ;
    ifstream trainingfile (argv[1]) ;
    if (trainingfile.is_open()){
        while ( trainingfile.good() ){
            getline (trainingfile,line) ;
            // delete text delimiter
            string format = "" ;
            boost::regex re("\"");
            if(boost::regex_search(line,re)){
                  line = boost::regex_replace(line,re,format,boost::format_default) ;
            }
            // split readline here :
            typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
            char_separator<char> sep(";");
            tokenizer tokens(line, sep);
            int i = 0 ;
            for(tokenizer::iterator tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter){
                if(i == MOLIDPOS){
                    TR_molid.push_back(*tok_iter) ;
                }else if(i == DESCRIPTORPOS){
                    vector<double> descrip_tmp  ;
                    string descriptors = *tok_iter ;
                    typedef boost::tokenizer<boost::char_separator<char> > tokenizer2 ;
                    char_separator<char> nosep(",")    ;
                    tokenizer2 tokens_tmp(descriptors, nosep);
                    for(tokenizer2::iterator tok_iter_tmp = tokens_tmp.begin(); tok_iter_tmp != tokens_tmp.end(); ++tok_iter_tmp){
                        int count = str2int(*tok_iter_tmp) ;
                        descrip_tmp.push_back(count) ;
                    }
                    TR_descriptors.push_back(descrip_tmp) ;
                    descrip_tmp.clear() ;
                }
                i++ ;
            }
        }
    }
    else{
        cout <<"ERROR opening file: " << argv[1] << endl ;
        exit (1);
    }

/*
 2   Read Test Set file
 */
    vector<string> TS_molid                ; // test set mol identifier
    vector<vector<double> > TS_descriptors ; // test set descriptors
    ifstream testfile (argv[2]) ;
    if (testfile.is_open()){
        while ( testfile.good() ){
            getline (testfile,line) ;
            // delete text delimiter
            string format = "" ;
            boost::regex re("\"");
            if(boost::regex_search(line,re)){
                  line = boost::regex_replace(line,re,format,boost::format_default) ;
            }
            // split readline here :
            typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
            char_separator<char> sep(";");
            tokenizer tokens(line, sep);
            int i = 0 ;
            for(tokenizer::iterator tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter){
                if(i == MOLIDPOS){
                    TS_molid.push_back(*tok_iter) ;
                }else if(i == DESCRIPTORPOS){
                    vector<double> descrip_tmp  ;
                    string descriptors = *tok_iter ;
                    typedef boost::tokenizer<boost::char_separator<char> > tokenizer2 ;
                    char_separator<char> nosep(",")    ;
                    tokenizer2 tokens_tmp(descriptors, nosep);
                    for(tokenizer2::iterator tok_iter_tmp = tokens_tmp.begin(); tok_iter_tmp != tokens_tmp.end(); ++tok_iter_tmp){
                        int count = str2int(*tok_iter_tmp) ;
                        descrip_tmp.push_back(count) ;
                    }
                    TS_descriptors.push_back(descrip_tmp) ;
                    descrip_tmp.clear() ;
                }
                i++ ;
            }
        }
    }
    else{
        cout <<"ERROR opening file: " << argv[2] << endl ;
        exit (1);
    }

/*
 3   compute <Di> for training moleules :
 */
    int nbrTRmol = TR_molid.size() ;
    vector<double> Di ;
    double Di_average = 0 ;
	for(int i=0 ; i < TR_descriptors.size() ; i++){
        vector<double> TR_tc_val ;
        for(int j=0 ; j < TR_descriptors.size() ; j++){
            double tc = tanimoto(TR_descriptors[i],TR_descriptors[j]);
            TR_tc_val.push_back(tc) ;
        }
        // get TR[i] nearest neighbors : do it by sorting TR_tc_val by descending order
        vector<string> TR_molid_copy = TR_molid ;
        int gap  = TR_molid.size() ;
        int swapped = 1 ;
        while(gap > 1 || swapped == 1){
            if(gap > 1){
                gap /= SHRINK ;
            }
            int z=0 ; swapped=0 ;
            while(z+gap < nbrTRmol){
                if(TR_tc_val[z] < TR_tc_val[z+gap]){
                    double tmp       = TR_tc_val[z]     ; string tmp2 = TR_molid_copy[z] ;
                    TR_tc_val[z]     = TR_tc_val[z+gap] ; TR_molid_copy[z] = TR_molid_copy[z+gap] ;
                    TR_tc_val[z+gap] = tmp              ; TR_molid_copy[z+gap] = tmp2 ;
                    swapped=1 ;
                }
                z += 1 ;
            }
        }
		// Calculate Di :
		double D = 0 ;
		for(int k=1 ; k<=K ; k++){ // don't begin from 0 because 0 is the same molecule as mol[i]
            D+=(1-TR_tc_val[k]) ;
            //view kNN for each training molecule
            //cout << k <<"\t"<<TR_molid_copy[k]<<"\t" <<TR_tc_val[k] <<endl ;
		}
 		D/=K ;
        Di.push_back(D) ;
        Di_average += D ;
    }
    Di_average /= nbrTRmol ;
    // calculate std dev for Di : sigma
    double sigma = 0 ;
    for(int i=0 ; i<nbrTRmol ; i++){
        sigma+=pow((Di[i]-Di_average),2) ;
    }
    double tmp = (double)1/nbrTRmol   ;
    sigma = sqrt(tmp*sigma) ;
    //cout << Di_average << "\t"<< sigma << endl ;

/*
 4   get k nearest neighbor for test set molecule :
 */
	for(int i=0 ; i < TS_descriptors.size() ; i++){
        vector<double> TS_tc_val ;
        for(int j=0 ; j < TR_descriptors.size() ; j++){
            double tc = tanimoto(TS_descriptors[i],TR_descriptors[j]);
            TS_tc_val.push_back(tc) ;
        }
        // get TR[i] nearest neighbors : do it by sorting TR_tc_val by descending order
        vector<string> TR_molid_copy = TR_molid ;
        int gap      = TR_molid.size() ;
        int nbrTSmol = TR_molid.size() ;
        int swapped = 1 ;
        while(gap > 1 || swapped == 1){
            if(gap > 1){
                gap /= SHRINK ;
            }
            int z=0 ; swapped=0 ;
            while(z+gap < nbrTSmol){
                if(TS_tc_val[z] < TS_tc_val[z+gap]){
                    double tmp       = TS_tc_val[z]     ;
                    TS_tc_val[z]     = TS_tc_val[z+gap] ;
                    TS_tc_val[z+gap] = tmp ;
                    swapped=1 ;
                }
                z += 1 ;
            }
        }
        double D = 0 ;
        for(int k=0 ; k<K ; k++){
            D+=(1-TS_tc_val[k]) ;
        }
        D/=K ;
        // Test :
        double rule = Di_average+(Z*sigma)  ;
        if(D <= rule){
            //cout << TS_molid[i] << "\tTrue" << endl ;
            cout << TS_molid[i] << "\tTrue\t" << D << "\t" << rule << endl ;
        }else{
            //cout << TS_molid[i] << "\tFalse" << endl ;
            cout << TS_molid[i] << "\tFalse\t" << D << "\t" << rule << endl ;
        }
	}

    return 0;
}



void sethelp(char *pname){
	printf("Usage\n\t%s TrainingSet.csv TestSet.csv k z\n",pname);
	printf("Description\n\tkNN Applicability Domain for fingerprints using tanimoto similarity\n");
	printf("\tFurther details in doi: 10.1021/jm010488u\n");
	printf("Author\n\tJamel Meslamani - Structural Chemogenomics Group - 2010\n") ;
}

int str2int(string strConvert) {
  int intReturn = atoi(strConvert.c_str());
  return(intReturn);
}

double tanimoto(vector<double> fp1,vector<double> fp2){
    double xij	= 0 ; double xii = 0 ; double xjj = 0 ;
    for(int i=0 ; i<fp1.size() ; i++){
        xij += (fp1[i] * fp2[i]) ;
        xii += (fp1[i] * fp1[i]) ;
        xjj += (fp2[i] * fp2[i]) ;
    }
    double tc = xij/(xii+xjj-xij) ;
	return(tc);
}
