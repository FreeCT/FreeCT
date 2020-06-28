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

#include <vector>

struct iterative_params{
    size_t Nx;
    size_t Ny;
    size_t Nz;
    double alpha;
    double beta;
    double lambda;
    double delta;
    int num_neighbors;
    double * weights;
    double weights_scale;
};

void   initialize_2d_weights(struct iterative_params * ip);
double quadratic(int curr_idx,struct iterative_params * ip      , double * recon_volume);
double edge_preserving(int curr_idx,struct iterative_params * ip, double * recon_volume);

