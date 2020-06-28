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

#pragma once

#include <iostream>
#include <iomanip>

inline void init_spinner(){
    std::cout << " ";    
}

inline void update_spinner(size_t i,size_t max){

    if (i!=0)
        std::cerr << "\b\b\b\b";            

    // Print a little spinner thing to show processing of long task:
    switch (i%4){
    case 0:{std::cerr << "\b|" ;break;}
    case 1:{std::cerr << "\b/" ;break;}
    case 2:{std::cerr << "\b-" ;break;}
    case 3:{std::cerr << "\b\\";break;}
    }

    int percentage=(i+1)*100/max;
    std::cerr << " " << std::setfill(' ') << std::setw(2) << percentage << "%";

}

inline void destroy_spinner(){
    std::cerr << "\b\b\b\b\b";
    std::cout << "\b \b";
}
