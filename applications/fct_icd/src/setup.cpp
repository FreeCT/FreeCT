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
#include <fstream>
#include <cmath>
#include <cstring>

//#include <fct/fct_read.h>
#include <fct/FreeCTRead.h>
//#include <ctbb/ctbb_read.h>

#include <stdio.h> // Bad form and hopefully we'll remove sometime. Needed for Compatibility with raw data reader libraries
#include "recon_structs.h"
#include "parse_config.h"
#include "setup.h"


#define PI 3.141592653589793238

#define string_check(TAG_NAME) if(rp.TAG_NAME.compare("")==0){                               \
        std::cout << #TAG_NAME" not properly set in configuration. Check parameter file.\n"; \
            exit_flag=true;}

#define double_check(TAG_NAME) if (rp.TAG_NAME==0.0){ \
        std::cout << #TAG_NAME" not properly set in configuration. Check parameter file.\n"; \
            exit_flag=true;}

#define size_t_check(TAG_NAME) if (rp.TAG_NAME==(size_t)0){             \
        std::cout << #TAG_NAME" not properly set in configuration. Check parameter file.\n"; \
            exit_flag=true;}

#define check_zero_allowed(TAG_NAME) if (rp.TAG_NAME==-1L){     \
        std::cout << #TAG_NAME" not properly set in configuration. Check parameter file.\n"; \
            exit_flag=true;}

#define check_with_default(TAG_NAME,DEFAULT_VALUE) if (rp.TAG_NAME==0.0){    \
    std::cout << #TAG_NAME" was not explicitly set in parameter file. Using default value: "#DEFAULT_VALUE << std::endl;\
    rp.TAG_NAME=DEFAULT_VALUE;}

#define check_string_with_default(TAG_NAME,DEFAULT_VALUE) if (strcmp(rp.TAG_NAME.c_str(),"")==0){ \
    std::cout << #TAG_NAME" was not explicitly set in parameter file. Using default value: "#DEFAULT_VALUE << std::endl;\
    rp.TAG_NAME=DEFAULT_VALUE;}

inline bool exists(const std::string& name){
    std::ifstream f(name.c_str());
    return f.good();
}

void allocate_recon_arrays(struct recon_params * rp,struct ct_data * data);
void     free_recon_arrays(struct recon_params * rp,struct ct_data * data);

struct recon_params configure_recon_params(std::string filename){
    // This function will configure and do first pass validation on
    // all of the reconstruction parameters necessary to accomplish a
    // reconstruction.  This ensures that we have the minimum required
    // information to process our data, however it should be noted
    // that this function does NOT guarantee accuracy of the data
    // passed in the configuration file and nothing can other than
    // careful user validation prior to running the program.

    // All reconstruction parameters are set and checked in this
    // function except for the following:
    //        rp.z_first_view
    //        rp.z_start
    //        rp.z_end
    //        rp.num_voxels_z    
    // These are set after the read of our raw data in the function
    // "load_raw_data"
    
    /*--- Read and parse the passed configuration file ---*/
    struct recon_params rp={};

    // For parameters that allow "0" values, set their default to be -1 and check against this value
    rp.FileType=-1;
    rp.FileSubType=-1;
    rp.RawOffset=-1;
    rp.Zffs=-1;
    rp.Phiffs=-1;
    rp.wfbp_initialize=-1;
    
    if (!exists(filename)){
        std::cout << "Parameter file not found. Exiting" << std::endl;
        exit(1);
    }
    
    parse_config(filename,&rp);

    /*--- Perform sanity checks and/or set defaults. ---*/
    // This step ensures that we have the minimum required parameter set for a
    // reconstruction
    bool exit_flag=false;

    // Required paths
    string_check(sinogram_path);
    string_check(output_dir);
    string_check(output_file);

    // Required Scanner Geometry
    double_check(acquisition_fov);
    size_t_check(num_views_per_turn_without_ffs);
    double_check(focal_spot_radius);
    double_check(source_detector_distance);
    double_check(anode_angle);
    size_t_check(n_channels);
    size_t_check(Nrows);
    double_check(axial_detector_spacing);
    double_check(axial_focal_spot_shift);
    double_check(center_channel_non_ffs);    
    double_check(center_row);
    double_check(transaxial_detector_spacing);
    double_check(transaxial_focal_spot_shift);
    double_check(fan_angle_increment);

    // Required Scan Specifics
    double_check(table_feed_per_rotation);

    // Recon Geometry
    check_with_default(recon_fov,25.0);
    check_with_default(nx,512);
    check_with_default(ny,512);
    double_check(z_start);
    double_check(z_end);
    check_with_default(slice_thickness,1.0);

    // Required Iterative parameters
    check_string_with_default(penalty,"quadratic");
    check_with_default(lambda,0.3);
    check_with_default(delta,0.005);
    size_t_check(num_iterations);

    // John's Addition's
    check_zero_allowed(FileType);
    check_zero_allowed(FileSubType);
    check_zero_allowed(RawOffset);
    size_t_check(Readings);
    check_zero_allowed(Zffs);
    check_zero_allowed(Phiffs);
    check_zero_allowed(wfbp_initialize);
    double_check(CollSlicewidth);

    // Deprecated
    // ========================================
    // double_check(fov_radius);
    // double_check(dx);
    // double_check(dy);
    // double_check(dz);
    // check_with_default(center_voxel_x,0.0);
    // check_with_default(center_voxel_y,0.0);

    // Non-macro-able other checks we need to perform
    // Penalty function
    if ((strcmp(rp.penalty.c_str(),"quadratic")!=0)&&
        (strcmp(rp.penalty.c_str(),"edge-preserving")!=0)){
        std::cout << "Requesting penalty function, \"" << rp.penalty << "\", is invalid." << std::endl;
        std::cout << "Valid choices are \"quadratic\" or \"edge-preserving\" (without quotes)" << std::endl;
        exit_flag=true;
    }

    if (exit_flag){
        std::cout << "Missing some required parameters.  Please review the parameter file. Exiting." << std::endl;
        exit(1);
    }

    /*--- Derive any remaining parameter values that we'll use ---*/
    // Mapping like this is bad form, but a holdover for compatibility w/ Stefano's code
    rp.center_voxel_x = (rp.nx-1.0)/2.0;
    rp.center_voxel_y = (rp.ny-1.0)/2.0;
    rp.dx = rp.recon_fov/(double)rp.nx;
    rp.dy = rp.recon_fov/(double)rp.ny;
    rp.dz = rp.slice_thickness;
    rp.fov_radius = rp.recon_fov/2.0;

    rp.voxel_size_x = rp.dx;
    rp.voxel_size_y = rp.dy;
    rp.voxel_size_z = rp.dz;
    rp.num_voxels_x = rp.nx;
    rp.num_voxels_y = rp.ny;
    rp.num_voxels_z = rp.nz;

    std::cout << "voxel_size_z (prior): " << rp.voxel_size_z << std::endl;

    // Note: Num_views_per_turn ACCOUNTS for flying focal spot.
    //       This will automatically adjust tube_angle_increment
    //       and num_views_for_system_matrix as needed.
    int n_ffs=(size_t)pow(2.0,(double)rp.Zffs)*(size_t)pow(2.0,(double)rp.Phiffs);
    rp.num_views_per_turn = rp.num_views_per_turn_without_ffs*(size_t)pow(2.0,(double)rp.Zffs)*(size_t)pow(2.0,(double)rp.Phiffs);
    rp.tube_angle_increment = 2.0*PI/rp.num_views_per_turn;
    rp.tube_z_increment   = rp.table_feed_per_rotation/rp.num_views_per_turn;
    //rp.voxel_size_z       = round(rp.voxel_size_z/rp.tube_z_increment)*rp.tube_z_increment; // Snap requested voxel size/slice thickness to one acceptable for the rotating slices algorithm
    rp.voxel_size_z       = round(rp.voxel_size_z/(n_ffs*rp.tube_z_increment))*(n_ffs*rp.tube_z_increment); // Snap requested voxel size/slice thickness to one acceptable for the rotating slices algorithm
    rp.views_per_slice    = (size_t)(rp.voxel_size_z/rp.tube_z_increment);

    std::cout << "voxel_size_z (after): " << rp.voxel_size_z << std::endl;
    
    rp.beam_width = (double)rp.Nrows*rp.axial_detector_spacing;

    rp.beam_width_at_roi_edge = rp.beam_width * (rp.focal_spot_radius+rp.fov_radius)/rp.source_detector_distance;
    std::cout << "Beam width: " << rp.beam_width_at_roi_edge << std::endl;
    
    rp.num_views_for_system_matrix = (size_t)(2.0*PI/(rp.table_feed_per_rotation*rp.tube_angle_increment)*rp.beam_width_at_roi_edge);
    // guarantee that we have an even number
    if (rp.num_views_for_system_matrix%2==1)
        rp.num_views_for_system_matrix+=1;
    // Guarantee divisible by 2*n_ffs
    while ((rp.num_views_for_system_matrix/2)%n_ffs > 0)
        rp.num_views_for_system_matrix=rp.num_views_for_system_matrix+1;

    rp.Nrows_projection=rp.Nrows/(size_t)pow(2.0,rp.Zffs);

    rp.center_voxel_z = 0.0;

    // Derived string parameters
    if (rp.matrix_path.compare("")==0)
        rp.matrix_path=rp.output_dir+"/matrix.bin";

    // We don't validate this at this stage (i.e. we don't know
    // whether or not z_end and z_start are actually in our dataset,
    // but we need an estimate of how many slices so that we can
    // allocate the reconstructed volume
    rp.num_voxels_z = (size_t)(fabs(rp.z_end-rp.z_start)/rp.voxel_size_z);
    std::cout << "Nvoxelsz= " << rp.num_voxels_z << std::endl;
    std::cout << "z end= " << rp.z_end << std::endl;
    std::cout << "z start= " << rp.z_start << std::endl;        

    return rp;
}

void load_raw(struct recon_params *  rp,struct ct_data * data){
    
    /*---  Ensure that the raw data file exists ---*/
    std::cout << "Checking for raw file: " << rp->sinogram_path << std::endl;
    if (!exists(rp->sinogram_path)){
        std::cout << "Raw data file not found. Exiting" << std::endl;
        exit(1);
    }
    std::cout << "Found!" << std::endl;

    std::cout << "Loading projection data..." << std::endl;
    
    allocate_recon_arrays(rp,data);

    /*--- Get tube angles and table positions ---*/
    // Bad form to use C version of file handling.  We'll port away at
    // some point, but current raw data reader libraries are compiled
    // to use this
    FILE * raw_file;
    raw_file=fopen(rp->sinogram_path.c_str(),"r");

    float * tmp_frame=new float [rp->Nrows_projection*rp->n_channels];

    int direction;
    
    switch (rp->FileType){
    case 0:{; // Binary file

	    for (int i=0;i<rp->Readings;i++){
                rp->first_view_angle=0.0;
                rp->table_direction=-1;
                
		data->tube_angles[i]=fmod(((360.0f/rp->num_views_per_turn)*i+rp->first_view_angle),360.0f);
		if (rp->table_direction==-1)
		    data->table_positions[i]=((float)rp->Readings/(float)rp->num_views_per_turn)
                        *rp->table_feed_per_rotation-(float)i*rp->table_feed_per_rotation/(float)rp->num_views_per_turn;
		else if (rp->table_direction==1)
		    data->table_positions[i]=0.0f+(float)i*rp->table_feed_per_rotation/(float)rp->num_views_per_turn;
		else 
		    data->table_positions[i]=0.0f+(float)i*rp->table_feed_per_rotation/(float)rp->num_views_per_turn;

                float * tmp_frame=(float*)malloc(rp->n_channels*rp->Nrows_projection*sizeof(float));

                ReadBinaryFrame(raw_file,i,rp->n_channels,rp->Nrows_projection, tmp_frame,rp->RawOffset);

                //Transpose all of our data into the frame
                for (int ii = 0; ii < rp->Nrows_projection; ++ii){
                    for (int jj = 0; jj < rp->n_channels; ++jj){
                        int idx_in=ii+jj*rp->Nrows_projection;
                        int idx_out=jj+ii*rp->n_channels;

                        data->raw[i*rp->n_channels*rp->Nrows_projection+idx_out]=tmp_frame[idx_in];
                    }
                }
	    }

            direction=rp->table_direction;
            
	    break;}
    case 1:{; //DefinitionAS Raw
            for (int i=0;i<rp->Readings;i++){
                // Extract tube angles and table positions
        	data->tube_angles[i]=(double)ReadPTRTubeAngle(raw_file,i,rp->n_channels,rp->Nrows_projection);
        	data->table_positions[i]=((double)ReadPTRTablePosition(raw_file,i,rp->n_channels,rp->Nrows_projection))/1000.0;
                ReadPTRFrame(raw_file,i,rp->n_channels,rp->Nrows_projection,&data->raw[i*rp->n_channels*rp->Nrows_projection]);
            }
     	
            // <0 is decreasing table position >0 is increasing
            direction=(data->table_positions[100]-data->table_positions[0])/fabs(data->table_positions[100]-data->table_positions[0]);
         
            break;}
    case 4:{; //IMA (can wrap any of the above (except binary)
            int raw_data_subtype=rp->FileSubType; // Determine if we're looking for PTR or CTD


      
            for (int i=0;i<rp->Readings;i++){
                data->tube_angles[i]=(double)ReadIMATubeAngle(raw_file,i,rp->n_channels,rp->Nrows_projection,raw_data_subtype,rp->RawOffset);
                data->table_positions[i]=((double)ReadIMATablePosition(raw_file,i,rp->n_channels,rp->Nrows_projection,raw_data_subtype,rp->RawOffset))/1000.0;
                ReadIMAFrame(raw_file,i,rp->n_channels,rp->Nrows_projection,&data->raw[i*rp->n_channels*rp->Nrows_projection],raw_data_subtype,rp->RawOffset);
            }

            // <0 is decreasing table position >0 is increasing
            direction=(data->table_positions[100]-data->table_positions[0])/fabs(data->table_positions[100]-data->table_positions[0]);
          
            break;}        
    }    
    fclose(raw_file);

    // Transpose all frames to match Stefano's implementation
    for (int i=0; i<rp->Readings; i++){
        std::memcpy(tmp_frame,&data->raw[i*rp->n_channels*rp->Nrows_projection],rp->n_channels*rp->Nrows_projection*sizeof(float));
        for (int j=0; j<rp->n_channels; j++){
            for (int k=0; k<rp->Nrows_projection; k++){
                size_t idx_in  = j+k*rp->n_channels;
                size_t idx_out = k+j*rp->Nrows_projection+i*rp->n_channels*rp->Nrows_projection;
                data->raw[idx_out]=tmp_frame[idx_in];
            }            
        }
    }
    
    delete[] tmp_frame;

    //TODO: need to parse table direction parameter in input file for use with non-Siemens data    
    rp->table_direction=direction;

    // Set z_first_view, the reset z_start, z_end, num_voxels_z
    rp->z_first_view = data->table_positions[0]/10.0; // convert mm to CM
    // Clean up the table positions (to get rid of the bad ones on the end)
    for (int i=0; i<rp->Readings; i++){
        data->table_positions[i]=rp->z_first_view+(double)i*(double)direction*rp->tube_z_increment;
    }

    // Check user requests against what we can actually reconstruct
    bool range_is_good=false;
    double allowed_start=rp->z_first_view+direction*ceil(0.5*rp->beam_width_at_roi_edge/rp->voxel_size_z)*rp->voxel_size_z;
    double allowed_end=data->table_positions[rp->Readings-1]-direction*ceil(0.5*rp->beam_width_at_roi_edge/rp->voxel_size_z)*rp->voxel_size_z;

    if (allowed_start>allowed_end){
        double tmp=allowed_start;
        allowed_start=allowed_end;
        allowed_end=tmp;
    }
    
    if ((rp->z_start>allowed_start) && (rp->z_start<allowed_end) &&
        (rp->z_end>allowed_start) && (rp->z_end<allowed_end))
        range_is_good=true;

    if (!range_is_good){
        std::cout << "Requested reconstruction range (" << rp->z_start << " : " << rp->z_end <<  ") is invalid" << std::endl;
        std::cout << "Allowed reconstruction range is: " << allowed_start << " : " << allowed_end << std::endl;
        exit(1);
    }

    rp->num_voxels_z = (size_t)(fabs(rp->z_end-rp->z_start)/rp->voxel_size_z);

    // Finally, go ahead an prepare our table of slice locations
    // and slice indices that we'll use for our slice rotations reconstructions
    data->slice_locations=new double [rp->num_voxels_z];
    data->slice_indices=new size_t [rp->num_voxels_z];

    int recon_direction=(rp->z_end-rp->z_start)/fabs(rp->z_end-rp->z_start);

    for (int i=0; i<rp->num_voxels_z; i++){
        // Location
        data->slice_locations[i]=rp->z_start+i*rp->voxel_size_z*recon_direction;

        // Index
        size_t central_idx=0;
        double curr_slice_location=data->slice_locations[i];
        int n_ffs=(size_t)pow(2.0,(double)rp->Zffs)*(size_t)pow(2.0,(double)rp->Phiffs);

        if (rp->table_direction>0){
            while (data->table_positions[central_idx]<curr_slice_location)
                central_idx+=n_ffs;                            
        }
        else{
            while (data->table_positions[central_idx]>curr_slice_location)
                central_idx+=n_ffs;                            
        }

        data->slice_indices[i]=central_idx;

        std::cout << "Slice " << i << ": " <<  data->slice_locations[i]  << " " << data->slice_indices[i] << std::endl;
    }

    std::cout << "voxel_size_z: "  << rp-> voxel_size_z  << std::endl;
    std::cout << "z_first_view: "  <<  rp->z_first_view  <<  std::endl;
    std::cout << "z_start: "       <<  rp->z_start       <<  std::endl;
    std::cout << "z_end: "         <<  rp->z_end         <<  std::endl;
    std::cout << "num_voxels_z: "  <<  rp->num_voxels_z  <<  std::endl;
    std::cout << "First_proj:  " <<  data->table_positions[0]              <<  std::endl;
    std::cout << "Last_proj:  "  <<  data->table_positions[rp->Readings-1] <<  std::endl;

    std::cout << "Done!" << std::endl;
}

void clean_up(struct recon_params *rp, struct ct_data *data){

    // Finalize the output reconstruction filepath (TODO: move the "else" statment up to "configure recon params")
    std::string output_filepath;
    if (rp->output_file!="")
        output_filepath = rp->output_dir+ "/" + rp->output_file;
    else
        output_filepath = rp->output_dir + "/" + "reconstruction_final.img";

    // Save our reconstructed volume to disk
    std::cout << "Writing final reconstruction to: " << output_filepath << std::endl;    
    std::ofstream recon_filepath(output_filepath,std::ios_base::binary);
    recon_filepath.write((char*)&data->recon_volume[0],rp->num_voxels_x*rp->num_voxels_y*rp->num_voxels_z*sizeof(float));
    recon_filepath.close();

    free_recon_arrays(rp,data);
}

void allocate_recon_arrays(struct recon_params * rp, struct ct_data * data){
    data->tube_angles     = new double [rp->Readings];
    data->table_positions = new double [rp->Readings];
    data->raw             = new float [rp->Readings*rp->Nrows_projection*rp->n_channels];
    data->recon_volume    = new float [rp->num_voxels_x*rp->num_voxels_y*rp->num_voxels_z](); // Initialize to zero
}

void free_recon_arrays(struct recon_params *rp, struct ct_data * data){
    delete[] data->tube_angles;
    delete[] data->table_positions;
    delete[] data->raw;
    delete[] data->recon_volume;
    delete[] data->slice_locations;
}
