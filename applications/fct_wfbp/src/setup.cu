/* FreeCT_wFBP is GPU and CPU CT reconstruction Software */
/* Copyright (C) 2015  John Hoffman */

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
/* jmhoffman@mednet.ucla.edu with "CTBANGBANG" in the subject line*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>

#include <setup.h>
#include <parse_config.h>
#include <fct/FreeCTRead.h>

#include <memory>

#define pi 3.1415926535897f
#define BLOCK_SLICES 32

void split_path_file(char**p, char**f, char *pf);
int array_search(float key,double * array,int numel_array,int search_type);
void remove_trailing_slash(char * str);

int configure_paths(struct recon_metadata *mr){
    
    /* --- Get working directory and User's home directory --- */
    struct passwd *pw=getpwuid(getuid());    
    //const char * home_dir=pw->pw_dir;
    //strcpy(mr->home_dir,home_dir);
    getcwd(mr->cwd,4096*sizeof(char));

    /* --- Get where the executable is running ---*/
    char full_exe_path[4096]={0};
    char * exe_path=(char*)calloc(4096,sizeof(char));
    char * exe_name=(char*)calloc(255,sizeof(char));
    readlink("/proc/self/exe",full_exe_path,4096);
    split_path_file(&exe_path,&exe_name,full_exe_path);
    strcpy(mr->install_dir,exe_path);
    mr->install_dir[strlen(mr->install_dir)-1]=0;

    /* --- Check for defined output path ---*/
    // if not defined, set to current working directory
    if (strcmp(mr->rp.output_dir,"")==0)
	strcpy(mr->rp.output_dir,mr->cwd);

    /* --- Check for output file name --- */
    // if not defined, set to rawdatafile.img
    if(strcmp(mr->rp.output_file,"")==0){
	char fullpath[4096+255]={0};
	sprintf(fullpath,"%s.img",mr->rp.raw_data_file);
	strcpy(mr->rp.output_file,fullpath);
    }

    // Cleanup directory strings
    //remove_trailing_slash(mr->home_dir);
    remove_trailing_slash(mr->install_dir);
    remove_trailing_slash(mr->cwd);    
    remove_trailing_slash(mr->rp.output_dir);
    remove_trailing_slash(mr->rp.raw_data_dir);    

    /* Check to make sure we can read the raw data file */
    char fullpath[4096+255]={0};
    FILE * fid;    
    memset(fullpath,0,4096+255);
    sprintf(fullpath,"%s/%s",mr->rp.raw_data_dir,mr->rp.raw_data_file);
    
    fid=fopen(fullpath,"r");
    if (fid==NULL){
	return 1;
    }
    else{
    	fclose(fid);
    }
    
    /* Check to make sure we can write to output file */
    memset(fullpath,0,4096+255);
    sprintf(fullpath,"%s/%s",mr->rp.output_dir,mr->rp.output_file);

    fid=fopen(fullpath,"w");
    if (fid==NULL){
	return 2;
    }
    else{
    	fclose(fid);
	// Theres a better way to do this... but for now this works
	remove(fullpath);	
    }
    
    return 0;
}

struct recon_params configure_recon_params(char * filename){
    struct recon_params prms;
    memset(&prms, 0,sizeof(prms));
    
    parse_config(filename,&prms);

    // Convert our table_dir_str to our table_dir integer
    if (strcmp(prms.table_dir_str,"")!=0){
	if (strcmp(prms.table_dir_str,"out")==0){
	    prms.table_dir=1;
	}
	else if (strcmp(prms.table_dir_str,"in")==0){
	    prms.table_dir=-1;
	}
	else{
	    printf("WARNING: TableDir parameter must be 'in' or 'out' (no quotes).  Defaulting to 'out'.\n");
	    prms.table_dir=1;
	}
    }

    // Perform some sanity checks to make sure that we have read in the "essentials"
    // Bail if critical values are zero
    int exit_flag=0;
    //if (prms.n_rows==0){
    //    printf("Nrows was not properly set in configuration.  Check parameter file.\n");
    //    exit_flag=1;
    //}
    //if (prms.coll_slicewidth==0){
    //    printf("CollSlicewidth was not properly set in configuration.  Check parameter file.\n");
    //    exit_flag=1;
    //}
    if (prms.slice_thickness==0){
	printf("SliceThickness was not properly set in configuration.  Check parameter file.\n");
	exit_flag=1;
    }
    //if (prms.pitch_value==0){
    //    printf("PitchValue was not properly set in configuration.  Check parameter file.\n");
    //    exit_flag=1;
    //}
    //if (prms.acq_fov==0){
    //    printf("AcqFOV was not properly set in configuration.  Check parameter file.\n");
    //    exit_flag=1;
    //}
    if (prms.recon_fov==0){
	printf("ReconFOV was not properly set in configuration.  Check parameter file.\n");
	exit_flag=1;
    }
    //if (prms.n_readings==0){
    //    printf("Readings was not properly set in configuration.  Check parameter file.\n");
    //    exit_flag=1;
    //}
    if (prms.nx==0){
	printf("Nx was not properly set in configuration.  Check parameter file.\n");
	exit_flag=1;
    }
    if (prms.ny==0){
	printf("Ny was not properly set in configuration.  Check parameter file.\n");
	exit_flag=1;
    }
    //if (prms.file_type==0&&prms.table_dir==0){
    //    printf("WARNING: 'TableDir' parameter unset.  Defaulting to 'out'.\n");
    //    prms.table_dir=1;
    //}
    if (exit_flag){
	exit(1);
    }
    
    return prms; 
} 

struct ct_geom configure_ct_geom(struct recon_metadata *mr){ 

    struct ct_geom cg;
    memset(&cg,0,sizeof(cg));

    // Runtime polymorphism to eventually support multiple raw data formats
    std::string raw_data_path = mr->rp.raw_data_dir;
    std::unique_ptr<fct::RawDataSet> ds = std::make_unique<fct::DicomDataSet>();
    ds->setPath(raw_data_path);
    ds->initialize();
    ds->readMetadata();

    // Physical geometry of the scanner (cannot change from scan to scan)
    // TCIA/Mayo clinic format does not automatically account for FFS in projections per rotation
    // Number of detector rows
    cg.anode_angle=7.0f*pi/180.0f; //!!!!!! How to get this for GE scanners?  Maybe we don't need it?
    cg.r_f             = ds->getDistSourceToIsocenter(); 
    cg.src_to_det      = ds->getDistSourceToDetector();
    cg.central_channel = ds->getDetectorCentralChannel();
    cg.n_rows          = ds->getDetectorRows();
    cg.n_channels      = ds->getDetectorChannels();
    cg.n_channels_oversampled = 2 * ds->getDetectorChannels();
    
    if (ds->getFlyingFocalSpotMode()=="FFSNONE"){
      mr->rp.phi_ffs = 0;
      mr->rp.z_ffs   = 0;
    }
    else if (ds->getFlyingFocalSpotMode()=="FFSXY"){
      mr->rp.phi_ffs = 1;
      mr->rp.z_ffs   = 0;
    }

    else if (ds->getFlyingFocalSpotMode()=="FFSZ"){
      mr->rp.phi_ffs = 0;
      mr->rp.z_ffs   = 1;
    }
    else if (ds->getFlyingFocalSpotMode()=="FFSXYZ"){
      mr->rp.phi_ffs = 1;
      mr->rp.z_ffs   = 1;
    }
    else {
      std::cout << "ERROR: Unsupported flying focal spot mode!" << std::endl;
      exit(1);
    }

    cg.n_proj_turn = ds->getProjectionsPerRotation();
    cg.n_proj_ffs = ds->getProjectionsPerRotation();
    std::cout << "Unclear if number of projections in Mayo Clinic format accounts for FFS!" << std::endl;
    //cg.n_proj_ffs  = cg.n_proj_turn*pow(2,mr->rp.phi_ffs)*pow(2,mr->rp.z_ffs);
    cg.n_rows_raw  = cg.n_rows;     //(unsigned int)(rp.n_rows/pow(2,rp.z_ffs));
    cg.n_rows      = cg.n_rows*pow(2,mr->rp.z_ffs);
    
    cg.fan_angle_increment = atan(ds->getDetectorTransverseSpacing()/ds->getDistSourceToDetector());
    mr->rp.coll_slicewidth = ds->getDistSourceToIsocenter()*(ds->getDetectorAxialSpacing()/ds->getDistSourceToDetector());

    // To accurately compute z_rot (table feed per rotation in mm)
    // we load 0->(n_rot) and read compute the travel
    ds->readProjection(0);
    ds->readProjection(cg.n_proj_ffs);

    cg.z_rot = ds->getTablePosition(cg.n_proj_ffs) - ds->getTablePosition(0);

    cg.add_projections     = (cg.fan_angle_increment*cg.n_channels/2)/(2.0f*pi/cg.n_proj_turn)+10; 
    cg.add_projections_ffs = cg.add_projections*pow(2,mr->rp.z_ffs)*pow(2,mr->rp.phi_ffs);

    //cg.theta_cone=2.0f*atan(7.5f*1.2f/cg.r_f);
    float detector_cone_offset = ((float)(cg.n_rows - 1))/2.0f;
    cg.theta_cone=2.0f*atan(detector_cone_offset * mr->rp.coll_slicewidth/cg.r_f);

    cg.acq_fov = 2.0f * cg.r_f*sin((float(cg.n_channels-1)/2.0) * ds->getDetectorTransverseSpacing() * (1.0f/cg.src_to_det));

    cg.table_direction = cg.z_rot/fabs(cg.z_rot);
    cg.z_rot = fabs(cg.z_rot);

    std::cout << "CT Geometry and Scan derived parameters: " << std::endl;
    std::cout << "===========================================" << std::endl;
    std::cout << "Num projections per turn:          "       << cg.n_proj_turn            << std::endl;
    std::cout << "Num projections per turn (ffs):    "       << cg.n_proj_ffs             << std::endl;
    std::cout << "Num detector channels:             "       << cg.n_channels             << std::endl;
    std::cout << "Num detector channels (ffs):       "       << cg.n_channels_oversampled << std::endl;
    std::cout << "Num detector rows:                 "       << cg.n_rows_raw             << std::endl;
    std::cout << "Num detector rows (ffs):           "       << cg.n_rows                 << std::endl;
    std::cout << "Radius src->isocenter (mm):        "       << cg.r_f                    << std::endl;
    std::cout << "Radius src->detector (mm):         "       << cg.src_to_det             << std::endl;
    std::cout << "Table feed per rotation (mm):      "       << cg.z_rot                  << std::endl;
    std::cout << "Theta cone (rad):                  "       << cg.theta_cone             << std::endl;
    std::cout << "Fan angle increment (rad):         "       << cg.fan_angle_increment    << std::endl;
    std::cout << "Anode angle (rad) *:               "       << cg.anode_angle            << std::endl;
    std::cout << "Central channel:                   "       << cg.central_channel        << std::endl;
    std::cout << "Acquisition FOV (mm):              "       << cg.acq_fov                << std::endl;
    std::cout << "Projection block buffer:           "       << cg.add_projections        << std::endl;
    std::cout << "Projection block buffer (ffs):     "       << cg.add_projections_ffs    << std::endl;
    std::cout << "Table direction *:                 "       << cg.table_direction        << std::endl;
    std::cout << "Collimated slicewidth @ isocenter: "       << mr->rp.coll_slicewidth    << std::endl;
    
    return cg;
}

void configure_reconstruction(struct recon_metadata *mr){
    /* --- Get tube angles and table positions --- */
    struct ct_geom cg=mr->cg;
    struct recon_params rp=mr->rp;

    // Allocate the memory
    mr->tube_angles=(float*)calloc(rp.n_readings,sizeof(float));
    mr->table_positions=(double*)calloc(rp.n_readings,sizeof(double));
    
    char fullpath[4096+255]={0};
    sprintf(fullpath,"%s/%s",rp.raw_data_dir,rp.raw_data_file);
    
    FILE * raw_file;
    raw_file=fopen(fullpath,"rb");
    if (raw_file==NULL){
	perror("Raw data file not found.");
	exit(1);	
    }
    
    //switch (rp.file_type){
    //case 0:{; // Binary file
    //        for (int i=0;i<rp.n_readings;i++){
    //    	mr->tube_angles[i]=fmod(((360.0f/cg.n_proj_ffs)*i+rp.tube_start_angle),360.0f);
    //    	if (cg.table_direction==-1)
    //    	    mr->table_positions[i]=((float)rp.n_readings/(float)cg.n_proj_ffs)*cg.z_rot-(float)i*cg.z_rot/(float)cg.n_proj_ffs;
    //    	else if (cg.table_direction==1)
    //    	    mr->table_positions[i]=0.0f+(float)i*cg.z_rot/(float)cg.n_proj_ffs;
    //    	else 
    //    	    mr->table_positions[i]=0.0f+(float)i*cg.z_rot/(float)cg.n_proj_ffs;
    //        }	
    //        break;}
    //case 1:{; //DefinitionAS Raw
    //        for (int i=0;i<rp.n_readings;i++){
    //    	mr->tube_angles[i]=ReadPTRTubeAngle(raw_file,i,cg.n_channels,cg.n_rows_raw);
    //    	mr->table_positions[i]=((double)ReadPTRTablePosition(raw_file,i,cg.n_channels,cg.n_rows_raw))/1000.0;		
    //        }
    //        
    //        // Clean up the table positions because they tend to
    //        // be wonky at the ends when read directly from the
    //        // raw data
    //    	
    //        // <0 is decreasing table position >0 is increasing
    //        int direction=(mr->table_positions[100]-mr->table_positions[0])/fabs(mr->table_positions[100]-mr->table_positions[0]);
    //        
    //        for (int i=1;i<rp.n_readings;i++){
    //    	mr->table_positions[i]=mr->table_positions[0]+(double)cg.z_rot*(((double)i)/(pow(2.0,rp.z_ffs)*pow(2.0,rp.phi_ffs)*(double)cg.n_proj_turn))*(double)direction;
    //        }
    //
    //        break;}
    //case 2:{; //CTD v1794 (Pre 2015 Sensation64)
    //        for (int i=0;i<rp.n_readings;i++){
    //    	mr->tube_angles[i]=ReadCTDv1794TubeAngle(raw_file,i,cg.n_channels,cg.n_rows_raw);
    //    	mr->table_positions[i]=(double)ReadCTDv1794TablePosition(raw_file,i,cg.n_channels,cg.n_rows_raw)/1000.0;
    //        }
    //        break;}
    //case 3:{; //CTD v2007 (Post 2015 Sensation64)
    //        for (int i=0;i<rp.n_readings;i++){
    //    	mr->tube_angles[i]=ReadCTDv2007TubeAngle(raw_file,i,cg.n_channels,cg.n_rows_raw);
    //    	mr->table_positions[i]=(double)ReadCTDv2007TablePosition(raw_file,i,cg.n_channels,cg.n_rows_raw)/1000.0;
    //        }
    //        break;}
    //case 4:{; //IMA (can wrap any of the above (except binary)
    //        int raw_data_subtype=mr->rp.file_subtype; // Determine if we're looking for PTR or CTD
    //    
    //        for (int i=0;i<rp.n_readings;i++){
    //    	mr->tube_angles[i]=ReadIMATubeAngle(raw_file,i,cg.n_channels,cg.n_rows_raw,raw_data_subtype,rp.raw_data_offset);
    //    	mr->table_positions[i]=((double)ReadIMATablePosition(raw_file,i,cg.n_channels,cg.n_rows_raw,raw_data_subtype,rp.raw_data_offset))/1000.0;
    //        }
    //
    //        // Clean up the table positions because they tend to
    //        // be wonky at the ends when read directly from the
    //        // raw data
    //
    //        // <0 is decreasing table position >0 is increasing
    //        int direction=(mr->table_positions[100]-mr->table_positions[0])/fabs(mr->table_positions[100]-mr->table_positions[0]);
    //        
    //        for (int i=1;i<rp.n_readings;i++){
    //    	mr->table_positions[i]=mr->table_positions[0]+(double)cg.z_rot*(((double)i)/(pow(2.0,rp.z_ffs)*pow(2.0,rp.phi_ffs)*(double)cg.n_proj_turn))*(double)direction;
    //        }
    //        
    //        break;}
    //case 5:{; //Force Raw
    //        for (int i=0;i<rp.n_readings;i++){
    //    	mr->tube_angles[i]=ReadForceTubeAngle(raw_file,i,cg.n_channels,cg.n_rows_raw);
    //    	mr->table_positions[i]=(double)ReadForceTablePosition(raw_file,i,cg.n_channels,cg.n_rows_raw)/1000.0;
    //        }
    //        break;}
    //case 6:{; //DICOM Raw
    //        for (int i=0;i<rp.n_readings;i++){
    //    	mr->tube_angles[i]=ReadDICOMTubeAngle(raw_file,i,cg.n_channels,cg.n_rows_raw);
    //    	mr->table_positions[i]=(double)ReadDICOMTablePosition(raw_file,i,cg.n_channels,cg.n_rows_raw)/1000.0;
    //        }
    //        break;}
    //}
    fclose(raw_file);

    /* --- Figure out how many and which projections to grab --- */

    int n_ffs=pow(2,rp.z_ffs)*pow(2,rp.phi_ffs);
    int n_slices_block=BLOCK_SLICES;

    int recon_direction=fabs(rp.end_pos-rp.start_pos)/(rp.end_pos-rp.start_pos);
    if (recon_direction!=1&&recon_direction!=-1) // user request one slice (end_pos==start_pos)
	recon_direction=1;

    // override end_pos if user has set the number of slices
    if (rp.n_slices!=0){
	rp.end_pos=rp.start_pos+(rp.n_slices-1)*rp.slice_thickness;
    }
    
    float recon_start_pos = rp.start_pos - recon_direction*rp.slice_thickness;
    float recon_end_pos   = rp.end_pos   + recon_direction*rp.slice_thickness;//rp.start_pos+recon_direction*(n_slices_recon-1)*rp.coll_slicewidth;

    int n_slices_requested=floor(fabs(recon_end_pos-recon_start_pos)/rp.coll_slicewidth)+1;//floor(fabs(rp.end_pos-rp.start_pos)/rp.coll_slicewidth)+1;
    int n_slices_recon=(n_slices_requested-1)+(n_slices_block-(n_slices_requested-1)%n_slices_block);

    recon_end_pos=recon_start_pos+recon_direction*(n_slices_recon-1)*rp.coll_slicewidth;
    
    int n_blocks=n_slices_recon/n_slices_block;

    //float recon_start_pos=rp.start_pos;
    //float recon_end_pos=rp.start_pos+recon_direction*(n_slices_recon-1)*rp.coll_slicewidth;
    int array_direction=fabs(mr->table_positions[100]-mr->table_positions[0])/(mr->table_positions[100]-mr->table_positions[0]);
    int idx_slice_start=array_search(recon_start_pos,mr->table_positions,rp.n_readings,array_direction);
    int idx_slice_end=array_search(recon_end_pos,mr->table_positions,rp.n_readings,array_direction);

    // Decide if the user has requested a valid range for reconstruction
    mr->ri.data_begin_pos = mr->table_positions[0];
    mr->ri.data_end_pos   = mr->table_positions[rp.n_readings-1];
    float projection_padding= cg.z_rot * (cg.n_proj_ffs/2+cg.add_projections_ffs+256)/cg.n_proj_ffs;
    float allowed_begin = mr->ri.data_begin_pos+array_direction*projection_padding;
    float allowed_end   = mr->ri.data_end_pos-array_direction*projection_padding;

    mr->ri.allowed_begin = allowed_begin;
    mr->ri.allowed_end   = allowed_end;

    // Check "testing" flag, write raw to disk if set
    if (mr->flags.testing){
	char fullpath[4096+255];
	strcpy(fullpath,mr->rp.output_dir);
	strcat(fullpath,"/table_positions.ct_test");
	FILE * outfile=fopen(fullpath,"w");
	fwrite(mr->table_positions,sizeof(double),rp.n_readings,outfile);
	fclose(outfile);

	strcpy(fullpath,mr->rp.output_dir);
	strcat(fullpath,"/tube_angles.ct_test");
	outfile=fopen(fullpath,"w");
	fwrite(mr->tube_angles,sizeof(float),rp.n_readings,outfile);
	fclose(outfile);
    }

    if (((rp.start_pos>allowed_begin)&&(rp.start_pos>allowed_end))||((rp.start_pos<allowed_begin)&&(rp.start_pos<allowed_end))){
	printf("Requested reconstruction is outside of allowed data range: %.2f to %.2f\n",allowed_begin,allowed_end);
	exit(1);
    }
    
    if (((rp.end_pos>allowed_begin)&&(rp.end_pos>allowed_end))||((rp.end_pos<allowed_begin)&&(rp.end_pos<allowed_end))){
	printf("Requested reconstruction is outside of allowed data range: %.2f to %.2f\n",allowed_begin,allowed_end);
	exit(1);
    }

    // We always pull projections in the order they occur in the raw
    // data.  If the end_pos comes before the start position in the
    // array, we use the end_pos as the "first" slice to pull
    // projections for.  This method will take into account the
    // ordering of projections with ascending or descending table
    // position, as well as any slice ordering the user requests.
    
    int idx_pull_start;
    int idx_pull_end;

    int pre_post_buffer=cg.n_proj_ffs/2;
    if (rp.z_ffs==1){
	pre_post_buffer=cg.n_proj_ffs/2;
    }
    
    if (idx_slice_start>idx_slice_end){
	idx_pull_start=idx_slice_end-pre_post_buffer-cg.add_projections_ffs;
	idx_pull_start=(idx_pull_start-1)+(n_ffs-(idx_pull_start-1)%n_ffs);
	idx_pull_end=idx_slice_start+pre_post_buffer+cg.add_projections_ffs;
	idx_pull_end=(idx_pull_end-1)+(n_ffs-(idx_pull_end-1)%n_ffs);
    }
    else{
	idx_pull_start=idx_slice_start-pre_post_buffer-cg.add_projections_ffs;
	idx_pull_start=(idx_pull_start-1)+(n_ffs-(idx_pull_start-1)%n_ffs);
	idx_pull_end=idx_slice_end+pre_post_buffer+cg.add_projections_ffs;
	idx_pull_end=(idx_pull_end-1)+(n_ffs-(idx_pull_end-1)%n_ffs);
    }

    idx_pull_end+=256;
   
    int n_proj_pull=idx_pull_end-idx_pull_start;
    
    // Ensure that we have a number of projections divisible by 128 (because GPU)
    n_proj_pull=(n_proj_pull-1)+(128-(n_proj_pull-1)%128);
    idx_pull_end=idx_pull_start+n_proj_pull;
    
    // copy this info into our recon metadata
    mr->cg.table_direction=array_direction;
    mr->rp.end_pos=rp.end_pos;
    mr->ri.n_ffs=n_ffs;
    mr->ri.n_slices_requested=n_slices_requested;
    mr->ri.n_slices_recon=n_slices_recon;
    mr->ri.n_slices_block=n_slices_block;
    mr->ri.n_blocks=n_blocks;
    mr->ri.idx_slice_start=idx_slice_start;
    mr->ri.idx_slice_end=idx_slice_end; 
    mr->ri.recon_start_pos=recon_start_pos;
    mr->ri.recon_end_pos=recon_end_pos;;
    mr->ri.idx_pull_start=idx_pull_start;
    mr->ri.idx_pull_end=idx_pull_end;
    mr->ri.n_proj_pull=n_proj_pull;

    /* --- Allocate our raw data array and our rebin array --- */
    mr->ctd.raw=(float*)calloc(cg.n_channels*cg.n_rows_raw*n_proj_pull,sizeof(float));
    mr->ctd.rebin=(float*)calloc(cg.n_channels_oversampled*cg.n_rows*(n_proj_pull-2*cg.add_projections_ffs)/n_ffs,sizeof(float));
    mr->ctd.image=(float*)calloc(rp.nx*rp.ny*n_slices_recon,sizeof(float));
}

void update_block_info(recon_metadata *mr){

    struct recon_info ri=mr->ri;
    struct recon_params rp=mr->rp;
    struct ct_geom cg=mr->cg;

    free(mr->ctd.raw);
    free(mr->ctd.rebin);
    
    /* --- Figure out how many and which projections to grab --- */
    int n_ffs=pow(2,rp.z_ffs)*pow(2,rp.phi_ffs);

    int recon_direction=fabs(rp.end_pos-rp.start_pos)/(rp.end_pos-rp.start_pos);
    if (recon_direction!=1&&recon_direction!=-1) // user requests one slice (end_pos==start_pos)
	recon_direction=1;
    
    float block_slice_start=ri.recon_start_pos+recon_direction*ri.cb.block_idx*rp.coll_slicewidth*(float)ri.n_slices_block;
    float block_slice_end=block_slice_start+(float)recon_direction*((float)ri.n_slices_block-1.0f)*rp.coll_slicewidth;
    int array_direction=fabs(mr->table_positions[100]-mr->table_positions[0])/(mr->table_positions[100]-mr->table_positions[0]);
    int idx_block_slice_start=array_search(block_slice_start,mr->table_positions,rp.n_readings,array_direction);
    int idx_block_slice_end=array_search(block_slice_end,mr->table_positions,rp.n_readings,array_direction);

    // We always pull projections in the order they occur in the raw
    // data.  If the end_pos comes before the start position in the
    // array, we use the end_pos as the "first" slice to pull
    // projections for.  This method will take into account the
    // ordering of projections with ascending or descending table
    // position, as well as any slice ordering the user requests.
    
    int idx_pull_start;
    int idx_pull_end;

    int pre_post_buffer=cg.n_proj_ffs/2;
    if (rp.z_ffs==1){
	pre_post_buffer=cg.n_proj_ffs/2;
    }

    if (idx_block_slice_start>idx_block_slice_end){
	idx_pull_start=idx_block_slice_end-pre_post_buffer-cg.add_projections_ffs;
	idx_pull_start=(idx_pull_start-1)+(n_ffs-(idx_pull_start-1)%n_ffs);
	idx_pull_end=idx_block_slice_start+pre_post_buffer+cg.add_projections_ffs;
	idx_pull_end=(idx_pull_end-1)+(n_ffs-(idx_pull_end-1)%n_ffs);
    }
    else{
	idx_pull_start=idx_block_slice_start-pre_post_buffer-cg.add_projections_ffs;
	idx_pull_start=(idx_pull_start-1)+(n_ffs-(idx_pull_start-1)%n_ffs);
	idx_pull_end=idx_block_slice_end+pre_post_buffer+cg.add_projections_ffs;
	idx_pull_end=(idx_pull_end-1)+(n_ffs-(idx_pull_end-1)%n_ffs);
    }

    idx_pull_end+=256;
   
    int n_proj_pull=idx_pull_end-idx_pull_start;

    // Ensure that we have a number of projections divisible by 128 (because GPU)
    n_proj_pull=(n_proj_pull-1)+(128-(n_proj_pull-1)%128);
    idx_pull_end=idx_pull_start+n_proj_pull;
    
    // copy this info into our recon metadata
    mr->ri.cb.block_slice_start=block_slice_start;
    mr->ri.cb.block_slice_end=block_slice_end;
    mr->ri.cb.idx_block_slice_start=idx_block_slice_start;
    mr->ri.cb.idx_block_slice_end=idx_block_slice_end; 

    mr->ri.idx_pull_start=idx_pull_start;
    mr->ri.idx_pull_end=idx_pull_end;
    mr->ri.n_proj_pull=n_proj_pull;

    mr->ri.cb.block_idx++;

    // Reallocate our raw and rebin arrays to account for changing n_proj_pull
    mr->ctd.raw=(float*)calloc(cg.n_channels*cg.n_rows_raw*n_proj_pull,sizeof(float));
    mr->ctd.rebin=(float*)calloc(cg.n_channels_oversampled*cg.n_rows*(n_proj_pull-2*cg.add_projections_ffs)/n_ffs,sizeof(float));
    
}

void extract_projections(struct recon_metadata * mr){

    float * frame_holder=(float*)calloc(mr->cg.n_channels*mr->cg.n_rows_raw,sizeof(float));

    FILE * raw_file;
    struct recon_params rp=mr->rp;
    struct ct_geom cg=mr->cg;
    char fullpath[4096+255]={0};
    sprintf(fullpath,"%s/%s",rp.raw_data_dir,rp.raw_data_file);
    raw_file=fopen(fullpath,"rb");
    
    //switch (mr->rp.file_type){
    //case 0:{ // binary
    //    for (int i=0;i<mr->ri.n_proj_pull;i++){
    //        ReadBinaryFrame(raw_file,mr->ri.idx_pull_start+i,cg.n_channels,cg.n_rows_raw,frame_holder,mr->rp.raw_data_offset);
    //        for (int j=0;j<cg.n_channels*cg.n_rows_raw;j++){
    //    	mr->ctd.raw[j+cg.n_channels*cg.n_rows_raw*i]=frame_holder[j];
    //        }
    //    }
    //    break;}
    //case 1:{ // DefinitionAS
    //    for (int i=0;i<mr->ri.n_proj_pull;i++){
    //        ReadPTRFrame(raw_file,mr->ri.idx_pull_start+i,cg.n_channels,cg.n_rows_raw,frame_holder);
    //        for (int j=0;j<cg.n_channels*cg.n_rows_raw;j++){
    //    	mr->ctd.raw[j+cg.n_channels*cg.n_rows_raw*i]=frame_holder[j];
    //        }
    //    }
    //    break;}
    //case 2:{ // CTD v1794 
    //    for (int i=0;i<mr->ri.n_proj_pull;i++){
    //        ReadCTDv1794Frame(raw_file,mr->ri.idx_pull_start+i,cg.n_channels,cg.n_rows_raw,frame_holder);
    //        for (int j=0;j<cg.n_channels*cg.n_rows_raw;j++){
    //    	mr->ctd.raw[j+cg.n_channels*cg.n_rows_raw*i]=frame_holder[j];
    //        }
    //    }
    //    break;}
    //case 3:{ // CTD v2007
    //    for (int i=0;i<mr->ri.n_proj_pull;i++){
    //        ReadCTDv2007Frame(raw_file,mr->ri.idx_pull_start+i,cg.n_channels,cg.n_rows_raw,frame_holder);
    //        for (int j=0;j<cg.n_channels*cg.n_rows_raw;j++){
    //    	mr->ctd.raw[j+cg.n_channels*cg.n_rows_raw*i]=frame_holder[j];
    //        }
    //    }
    //    break;}
    //case 4:{ // IMA (wraps either PTR or IMA)
    //    int raw_data_subtype=rp.file_subtype;
    //    for (int i=0;i<mr->ri.n_proj_pull;i++){
    //        ReadIMAFrame(raw_file,mr->ri.idx_pull_start+i,cg.n_channels,cg.n_rows_raw,frame_holder,raw_data_subtype,rp.raw_data_offset);
    //        for (int j=0;j<cg.n_channels*cg.n_rows_raw;j++){
    //    	mr->ctd.raw[j+cg.n_channels*cg.n_rows_raw*i]=frame_holder[j];
    //        }
    //    }
    //    break;}	
    //case 5:{ //Force Raw
    //    for (int i=0;i<mr->ri.n_proj_pull;i++){
    //        
    //        ReadForceFrame(raw_file,mr->ri.idx_pull_start+i,cg.n_channels,cg.n_rows_raw,frame_holder);
    //
    //        for (int j=0;j<cg.n_channels*cg.n_rows_raw;j++){
    //    	mr->ctd.raw[j+cg.n_channels*cg.n_rows_raw*i]=frame_holder[j];
    //        }
    //
    //    }
    //    break;}
    //case 6:{ //DICOM Raw
    //    for (int i=0;i<mr->ri.n_proj_pull;i++){
    //        ReadDICOMFrame(raw_file,mr->ri.idx_pull_start+i,cg.n_channels,cg.n_rows_raw,frame_holder);
    //        for (int j=0;j<cg.n_channels*cg.n_rows_raw;j++){
    //    	mr->ctd.raw[j+cg.n_channels*cg.n_rows_raw*i]=frame_holder[j];
    //        }
    //    }
    //    break;}
    //}

    // Check "testing" flag, write raw to disk if set
    if (mr->flags.testing){
	memset(fullpath,0,4096+255);
	strcpy(fullpath,mr->rp.output_dir);
	strcat(fullpath,"/raw.ct_test");
	FILE * outfile=fopen(fullpath,"w");
	fwrite(mr->ctd.raw,sizeof(float),cg.n_channels*cg.n_rows_raw*mr->ri.n_proj_pull,outfile);
	fclose(outfile);
    }
    
    fclose(raw_file);
    free(frame_holder);
}

void finish_and_cleanup(struct recon_metadata * mr){

    int n_slices_final=floor(fabs(mr->rp.end_pos-mr->rp.start_pos)/mr->rp.slice_thickness)+1;
    
    // Write the image data to disk
    char fullpath[4096+255]={0};
    sprintf(fullpath,"%s/%s",mr->rp.output_dir,mr->rp.output_file);
    FILE * outfile=fopen(fullpath,"w");
    fwrite(mr->ctd.final_image_stack,sizeof(float),mr->rp.nx*mr->rp.ny*n_slices_final,outfile);
    fclose(outfile);

    // Free all remaining allocations in metadata
    free(mr->ctd.rebin);
    free(mr->ctd.image);
    free(mr->ctd.raw);
    free(mr->ctd.final_image_stack);    
    free(mr->tube_angles);
    free(mr->table_positions);
}


void remove_trailing_slash(char * str){
    size_t len=strlen(str);
    if ((len>0)&&(str[len-1]=='/')){
	str[len-1]='\0';
    }
}

void split_path_file(char**p, char**f, char *pf) {
    char *slash = pf, *next;
    while ((next = strpbrk(slash + 1, "\\/"))) slash = next;
    if (pf != slash) slash++;
    *p = strndup(pf, slash - pf);
    *f = strdup(slash);
}


int array_search(float key,double * array,int numel_array,int search_type){
    int idx=0;

    switch (search_type){
    case -1:{// Array descending
	while (key<array[idx]&&idx<numel_array){
	    idx++;}
	break;}
    case 0:{// Find where we're equal
	while (key!=array[idx]&&idx<numel_array){
	    idx++;}
	break;}
    case 1:{// Array ascending
	while (key>array[idx]&&idx<numel_array){
	    idx++;}
	break;}
    }

    return idx;
}
