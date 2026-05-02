#pragma once

namespace qe {

int handle_cif_mode(int argc, char** argv, int s = 0);
int handle_dos_mode(int argc, char** argv, int s = 0);

int handle_band_pre_mode(int argc, char** argv, int s = 0);
int handle_band_post_mode(int argc, char** argv, int s = 0);
int handle_band_fat_mode(int argc, char** argv, int s = 0);

int handle_elastic_pre_mode(int argc, char** argv, int s = 0);
int handle_elastic_post_mode(int argc, char** argv, int s = 0);

int handle_charge_pre_mode(int argc, char** argv, int s = 0);
int handle_charge_post_mode(int argc, char** argv, int s = 0);

int handle_mag_mode(int argc, char** argv, int s = 0);

int handle_stm_pre_mode(int argc, char** argv, int s = 0);
int handle_stm_post_mode(int argc, char** argv, int s = 0);

int handle_bader_mode(int argc, char** argv, int s = 0);

}  // namespace qe
