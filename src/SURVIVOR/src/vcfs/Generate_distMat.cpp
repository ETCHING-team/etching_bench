/*
 * Generate_distMat.cpp
 *
 *  Created on: Jul 17, 2017
 *      Author: sedlazec
 */

#include "Generate_distMat.h"
bool is_file_exist(const char *fileName) {
	std::ifstream infile(fileName);
	return infile.good();
}

std::vector<std::string> parse_sample_names(std::string filename) {
	size_t buffer_size = 200000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;

	myfile.open(filename.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "VCF Parser: could not open file: " << filename.c_str() << std::endl;
		exit(0);
	}

	myfile.getline(buffer, buffer_size);

	std::vector<std::string> names;
	while (!myfile.eof()) {
		if (buffer[0] == '#' && buffer[1] != '#') {
			//parse line!
			int count = 0;
			std::string name = "";
			for (size_t i = 0; i < buffer_size && buffer[i] != '\n' && buffer[i] != '\0'; i++) {
				if (count > 8 && buffer[i] != '\t' && name.size()<64) {
					name += buffer[i];
				}
				if (buffer[i] == '\t') {
					if (!name.empty()) {
						names.push_back(name);
						name.clear();
					}
					count++;
				}
			}

			if (!name.empty()) {
				names.push_back(name);
				name.clear();
			}
		} else if (buffer[0] != '#') {
			break;
		}
		myfile.getline(buffer, buffer_size);
	}
	myfile.close();
	return names;
}
bool not_set(char * buffer){
	size_t i=0;
	while(buffer[i]!='\t' && buffer[i]!='\0'){
		if(strncmp("NaN",&buffer[i],3)==0){
			return true;
		}
		i++;
	}
	return false;
}
void update_mat(std::string filename, std::vector<std::vector<int> > &samples_mat, std::vector<std::string> sample_names) {
	//std::cerr<<"Parser was adapted for bug with two tabs!"<<std::endl;
	size_t buffer_size = 200000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;

	myfile.open(filename.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "VCF Parser: could not open file: " << filename.c_str() << std::endl;
		exit(0);
	}
	myfile.getline(buffer, buffer_size);
	int num=0;
	while (!myfile.eof()) {
		if (buffer[0] != '#') {
			int count = 0;
			int id = 0;
			//collect a vector ID that share
			std::vector<int> ids;
			for (size_t i = 0; i < buffer_size && buffer[i] != '\n' && buffer[i] != '\0'; i++) {
				if (count > 8 && buffer[i - 1] == '\t') {
					if (!not_set(&buffer[i]) && !(buffer[i] == '0' && buffer[i + 2] == '0')) {
						ids.push_back(id);
					}
					id++;

				}
				if (buffer[i] == '\t') {
					count++;
				}
			}
			//std::cout<<"ID: "<<id <<" maxx: "<<samples_mat.size()<<std::endl;
			//update mat:
			for (size_t i = 0; i < ids.size(); i++) {
				for (size_t j = 0; j < ids.size(); j++) { //check me!
					samples_mat[ids[i]][ids[j]]++;
				}
			}
			if(num%10000==0){
				std::cout<<" Processed SV: "<<num<<std::endl;
			}
			num++;
		}
		myfile.getline(buffer, buffer_size);
	}
}

void generate_dist_mat(std::string svs_vcf, std::string snp_vcf, std::string weighted_file, std::string output) {
	//determine # samples:
	std::vector<std::string> sample_names;
	if (is_file_exist(svs_vcf.c_str())) {
		sample_names = parse_sample_names(svs_vcf);
	} else if (is_file_exist(snp_vcf.c_str())) {
		sample_names = parse_sample_names(snp_vcf);
	} else {
		std::cerr << "We need at least one SNP/SVs file with all the samples!" << std::endl;
		exit(1);
	}
	std::cout << "We detected " << sample_names.size() << " Samples" << std::endl;

	//initialize:
	std::vector<std::vector<int> > samples_mat;
	std::vector<int> tmp;
	tmp.assign(sample_names.size(), 0);
	samples_mat.assign(sample_names.size(), tmp);

	std::cout<<"Finished int"<<std::endl;

	if (is_file_exist(svs_vcf.c_str())) {
		update_mat(svs_vcf, samples_mat, sample_names);
	}
	if (is_file_exist(snp_vcf.c_str())) {
		update_mat(snp_vcf, samples_mat, sample_names); //TODO test for SNP!
	}
	if (is_file_exist(weighted_file.c_str())) {
		//	update_mat_weights(weighted_file, samples_mat, sample_names);
	}
	std::cout<<"Printing:"<<std::endl;
	//print matrix:
	FILE * file = fopen(output.c_str(), "w");
	fprintf(file, "%i",(int)sample_names.size());
	fprintf(file, "%c", '\n');

	for (size_t i = 0; i != sample_names.size(); i++) {
		fprintf(file, "%s",sample_names[i].c_str());
		for (size_t j = 0; j <sample_names.size(); j++) {
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", samples_mat[i][j]);
		}
		fprintf(file, "%c", '\n');
	}
	fclose(file);
}
