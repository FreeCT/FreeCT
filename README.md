# FreeCT
FreeCT is open source reconstruction software for clinical 3rd generation, fan-beam CT

This repository is a temporary placeholder for what will eventually become the monolithic FreeCT package.  Until it is converted over, we will use it as a unified information gathering point about converting the project to run on the TCIA LDCT raw projection dataset (https://wiki.cancerimagingarchive.net/pages/viewpage.action?pageId=52758026).

If you,

 * are interested in contributing to accelerate this process,
 * have feature requests, or 
 * ideas for how to improve the project
 
please reach out to John at johnmarianhoffman@gmail.com.

All "timelines" provided below are best guess estimates and are subject to change.

## Near-term work to enable use on TCIA datasets:
Short-term goal is to get SOMETHING running for the LDCT dataset.  This primarily involves modifications to the reader library and "setup" components of the FreeCT_wFBP packages.  ICD will eventually be ported, however we do not have a timeline for that.

Timeline for work: 2 weeks - 1 month 

* Convert FreeCT_Reader to support TCIA DICOM format (**done**)
* Modify FreeCT_wFBP setup.cu: support reading scanner data directly from RawDataSet class (*in progress*)
* Modify FreeCT_wFBP setup.cu: New structure for YAML input files
* Any required modifications to the reconstruction process b/c of new input structures

(updated on 2020/06/28)

## Mid-term work to improve FreeCT usability
Larger, project-level restructuring to make downloading, building, and installing required FreeCT components easier.  This will necessitate a license change (likely going to Apache 2.0 license.)

Timeline for work: 2 - 6 months 

* Migrate all project components (Reader, WFBP, ICD) into a single repository (this repository) with a unified build process
* Unified configuration file structure/parsing framework
* Packaging of the project Ubuntu (and potentially other high-use OSes)

(updated on 2020/06/28)

## Long-term work
Timeline for work: 6 months-???

* Modernize GPU code for larger-memory cards
* Better decoupling of CUDA code from C++ code
* Python wrappers
* Interoperability with other reconstruction software packages

(updated on 2020/06/28)
