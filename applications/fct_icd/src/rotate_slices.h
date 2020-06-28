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

// Rotate slices from fixed cartesian grid into rotating coordinate frame
// e.g. Initial WFBP slices into ICD rotating frame
void rotate_slices_fixed2rotating(const struct recon_params * rp,struct ct_data * data);

// Return rotating slices to fixed cartesian grid
// e.g. Reconstructed ICD slices back for final save to disk
void rotate_slices_rotating2fixed(const struct recon_params * rp,struct ct_data * data);
