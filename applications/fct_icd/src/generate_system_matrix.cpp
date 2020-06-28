/* FreeCT_ICD is MBIR CT reconstruction Software */
/* Copyright (C) 2018  John Hoffman, Stefano Young, Frederic Noo */

/* This program is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU General Public License */
/* as published by the Free Software Foundation; either version 2 */
/* of the License, or (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with this program; if not, write to the Free Software */
/* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. */

/* Questions and comments should be directed to */
/* jmhoffman@mednet.ucla.edu with "FreeCT_ICD" in the subject line*/


#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdio>
#include "generate_system_matrix.h"
#include "recon_structs.h"

#include <vector>

//#include <Eigen/Dense>
//#include <Eigen/SparseCore>

#include "spinner.h"

#include "boost/numeric/ublas/vector_sparse.hpp" 
#include "boost/numeric/ublas/vector.hpp"
#include "boost/numeric/ublas/io.hpp"
namespace ublas = boost::numeric::ublas;

#define PI 3.14159265359

#define debug_disp(VARIABLE_NAME) std::cout << #VARIABLE_NAME": " << VARIABLE_NAME << std::endl;

struct pair{
    int index;
    float value;
};

void save_system_matrix_first_block(const struct recon_params * rp,std::vector<ublas::compressed_vector<float>> & system_matrix_block);
void       save_system_matrix_block(const struct recon_params * rp,std::vector<ublas::compressed_vector<float>> & system_matrix_block);

void generate_system_matrix(const struct recon_params * rp, struct ct_data * data){
    std::cout << "Generating system matrix (no flying focal spots)..." << std::endl;

    ublas::vector<double> source_position(3);
    ublas::vector<double> direction_cosine(3);
    ublas::vector<double> e_w(3);
    ublas::vector<double> e_u(3);
    ublas::vector<double> e_z(3);

    std::cout << rp->num_voxels_x << " : " << rp->num_voxels_y << std::endl;
    
    std::vector<ublas::compressed_vector<float>> system_matrix(rp->num_voxels_x*rp->num_voxels_y);
    for (int i = 0; i < system_matrix.size(); i++){
        system_matrix[i].resize((size_t)(rp->Readings*rp->n_channels*rp->Nrows_projection));
    }

    // Report initial information about system matrix
    std::cout << "System matrix size: " << system_matrix.size() << std::endl;
    std::cout << "Allocating system matrix..." << std::endl;

    // Attempts with BOOST
    if ((!rp->Zffs)&(!rp->Phiffs)){
        std::cout << "Matrix generation for SciFFSNone configuration."  << std::endl;

        init_spinner();

        for (int i = 0; i < rp->num_views_for_system_matrix; i++){
            
            update_spinner(i,rp->num_views_for_system_matrix);

            // tk this could be one place where we take the initial tube angle into account
            double tube_angle = (i - 0.5*rp->num_views_for_system_matrix)*rp->tube_angle_increment;
 
            //Define local coordinate system for the current view
            e_w(0) = cos(tube_angle);
            e_w(1) = sin(tube_angle);
            e_w(2) = 0.0;

            e_u(0) = -sin(tube_angle);
            e_u(1) = cos(tube_angle);
            e_u(2) = 0.0;

            e_z(0) = 0.0;
            e_z(1) = 0.0;
            e_z(2) = 1.0;

            source_position = rp->focal_spot_radius*e_w + rp->table_direction*((double)i - 0.5*(double)rp->num_views_for_system_matrix)*rp->tube_z_increment*e_z;

            for (int j = 0; j < rp->n_channels; j++){
                double transaxial_position = (j - rp->center_channel_non_ffs)*rp->transaxial_detector_spacing;
                double transaxial_angle = transaxial_position / rp->source_detector_distance;
                double cos_transaxial_angle = cos(transaxial_angle);
                double sin_transaxial_angle = sin(transaxial_angle); 

                for (int k = 0; k < rp->Nrows_projection; k++){
                    double axial_position = (k - rp->center_row)*rp->axial_detector_spacing;

                    direction_cosine = (-rp->source_detector_distance*cos_transaxial_angle*e_w  \
                                        -rp->source_detector_distance*sin_transaxial_angle*e_u  \
                                        -axial_position*e_z)/sqrt(axial_position*axial_position \
                                        +rp->source_detector_distance*rp->source_detector_distance);

                    int q = k + rp->Nrows_projection*j + rp->Nrows_projection*rp->n_channels*i;
                    //int q = (rp->Nrows_projection-1-k) + rp->Nrows_projection*j + rp->Nrows_projection*rp->n_channels*i;
                    
                    //Compute the contribution of the ray to the system matrix
                    double x, y,
                        x_hat, y_hat, z_hat,
                        i_hat, j_hat, k_hat;

                    double abs_alpha_x = fabs(direction_cosine(0));
                    double abs_alpha_y = fabs(direction_cosine(1));

                    if (abs_alpha_x > abs_alpha_y){

                        //std::cout << "x focused" << std::endl;
                        
                        for (int i = 0; i < rp->num_voxels_x; i++){
                            x = (i - rp->center_voxel_x)*rp->voxel_size_x;

                            y_hat = source_position(1) + direction_cosine(1) / direction_cosine(0)*(x - source_position(0));
                            j_hat = y_hat / rp->voxel_size_y + rp->center_voxel_y;
                            int j = (int)(floor(j_hat));

                            z_hat = source_position(2) + direction_cosine(2) / direction_cosine(0)*(x - source_position(0));
                            k_hat = z_hat / rp->voxel_size_z + rp->center_voxel_z;
                            int k = (int)(floor(k_hat)); // for negative values of k, casting to int gives the wrong answer so I use floor

                            if ((j >= 0) && (j <= (int)(rp->num_voxels_y - 2))){
                                double scale = rp->voxel_size_x / abs_alpha_x;
                                
                                if (k == -1){
                                    system_matrix[i + rp->num_voxels_x*j].push_back(q, (float)(scale*(k_hat - k)*(1 - (j_hat - j))));
                                    system_matrix[i + rp->num_voxels_x*(j + 1)].push_back(q, (float)(scale*(k_hat - k)*(j_hat - j)));                                    
                                }
                                if (k == 0){
                                    system_matrix[i + rp->num_voxels_x*j].push_back(q, (float)(scale*(1 - (k_hat - k))*(1 - (j_hat - j))));
                                    system_matrix[i + rp->num_voxels_x*(j + 1)].push_back(q, (float)(scale*(1 - (k_hat - k))*(j_hat - j)));
                                }
                            }
                        }
                    }
                    else{

                        //std::cout << "y focused" << std::endl;
                        
                        for (int j = 0; j < rp->num_voxels_y; j++){                            
                            y = (j - rp->center_voxel_y)*rp->voxel_size_y;

                            x_hat = source_position(0) + direction_cosine(0) / direction_cosine(1)*(y - source_position(1));
                            i_hat = x_hat / rp->voxel_size_x + rp->center_voxel_x;
                            int i = (int)(floor(i_hat));

                            z_hat = source_position(2) + direction_cosine(2) / direction_cosine(1)*(y - source_position(1));
                            k_hat = z_hat / rp->voxel_size_z + rp->center_voxel_z;
                            int k = (int)(floor(k_hat));

                            if ((i >= 0) && (i <= (int)(rp->num_voxels_x - 2))){
                                double scale = rp->voxel_size_y / abs_alpha_y;

                                if (k == -1){
                                    system_matrix[i +       rp->num_voxels_x*j].push_back(q, (float)(scale*(k_hat - k)*(1 - (i_hat - i))));
                                    system_matrix[(i + 1) + rp->num_voxels_x*j].push_back(q, (float)(scale*(k_hat - k)*(i_hat - i)));                                    
                                }
                                if (k == 0){
                                    system_matrix[i +       rp->num_voxels_x*j].push_back(q, (float)(scale*(1 - (k_hat - k))*(1 - (i_hat - i))));
                                    system_matrix[(i + 1) + rp->num_voxels_x*j].push_back(q, (float)(scale*(1 - (k_hat - k))*(i_hat - i)));
                                }
                            }
                        }
                    }
                }
            }

            // Save system matrix chunk to disk (should happen ~5 times throughout generation)
            if (i==300){
                save_system_matrix_first_block(rp,system_matrix);
            }
            else if (i%300==0 && i!=0){
                save_system_matrix_block(rp,system_matrix);
            }

        }

        // Save the final chunk to disk
        if (rp->num_views_for_system_matrix<300)
            save_system_matrix_first_block(rp,system_matrix);
        else            
            save_system_matrix_block(rp,system_matrix);
        
        // Clear our user feeback
        destroy_spinner();

        std::cout << "Done!" << std::endl;
    }
}


void save_system_matrix_first_block(const struct recon_params * rp,std::vector<ublas::compressed_vector<float>> & system_matrix_block){

    std::cout << "Saving first block to disk..." << std::endl;

    // Open the output file
    std::string matrix_filepath=rp->output_dir+"/matrix.bin";    
    std::ofstream output_file(matrix_filepath.c_str(),std::ios_base::binary);

    // Loop over all columns of our system matrix
    for (int i=0; i<rp->num_voxels_x*rp->num_voxels_y; i++){

        // Get the information about the column chunk in memory
        size_t num_nonzeros_memory=(size_t)system_matrix_block[i].filled();
        struct pair * col_memory=new struct pair[num_nonzeros_memory];
        
        size_t nonzero_idx=0;
        for (auto pos = system_matrix_block[i].begin(); pos != system_matrix_block[i].end(); ++pos){
            col_memory[nonzero_idx].index=(int)pos.index();
            col_memory[nonzero_idx].value=(*pos);
            nonzero_idx=nonzero_idx + 1;
        }

        // Flush the data to disk (more lines of code than needed for clarity);
        size_t num_nonzeros=num_nonzeros_memory;
        output_file.write((char*)&num_nonzeros,sizeof(num_nonzeros));
        output_file.write((char*)col_memory,num_nonzeros_memory*sizeof(struct pair));

        // Empty contents, and force reallocation of matrix column
        ublas::compressed_vector<float>(system_matrix_block[i].size(), 0).swap(system_matrix_block[i]);

        delete[] col_memory;
    }

    output_file.close();
}


void save_system_matrix_block(const struct recon_params * rp,std::vector<ublas::compressed_vector<float>> & system_matrix_block){

    std::cout << "Saving next block to disk..." << std::endl;

    // Filepaths
    std::string tmp_matrix_filepath=rp->output_dir + "/matrix_tmp.bin";
    std::string matrix_filepath=rp->output_dir+"/matrix.bin";

    // Open our input and output files
    std::ofstream   output_file(tmp_matrix_filepath.c_str(),std::ios_base::binary);
    std::ifstream existing_file(matrix_filepath.c_str(),std::ios_base::binary);

    for (int i=0; i<rp->num_voxels_x*rp->num_voxels_y; i++){

        // Read the already saved chunk of the current column from disk
        size_t num_nonzeros_existing;
        existing_file.read((char*)&num_nonzeros_existing,sizeof(num_nonzeros_existing));

        struct pair * col_disk = new struct pair[num_nonzeros_existing];
        existing_file.read((char*)col_disk,num_nonzeros_existing*sizeof(struct pair));

        // Get the information about the chunk in memory
        size_t num_nonzeros_memory=(size_t)system_matrix_block[i].filled();
        struct pair * col_memory=new struct pair[num_nonzeros_memory];
        
        size_t nonzero_idx=0;
        for (auto pos = system_matrix_block[i].begin(); pos != system_matrix_block[i].end(); ++pos){            
            col_memory[nonzero_idx].index=(int)pos.index();
            col_memory[nonzero_idx].value=(*pos);
            nonzero_idx=nonzero_idx + 1;
        }

        // Flush the data to disk (more lines of code than needed for clarity);
        size_t num_nonzeros=num_nonzeros_existing+num_nonzeros_memory;
        output_file.write((char*)&num_nonzeros,sizeof(num_nonzeros));
        output_file.write((char*)col_disk,num_nonzeros_existing*sizeof(struct pair));
        output_file.write((char*)col_memory,num_nonzeros_memory*sizeof(struct pair));

        // Empty contents, and force reallocation of matrix column
        ublas::compressed_vector<float>(system_matrix_block[i].size(), 0).swap(system_matrix_block[i]);

        delete[] col_memory;
        delete[] col_disk;
    }

    // Overwrite the existing matrix file with the updated temporary file
    int success=rename(tmp_matrix_filepath.c_str(),matrix_filepath.c_str()); // 0 return means success
    if (success!=0){
        std::cout << "There was an error joining matrix chunks. Exiting."  << std::endl;
        exit(1);
    }

    output_file.close();
    existing_file.close();
}
