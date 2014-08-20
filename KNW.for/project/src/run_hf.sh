Nmax=4
num_part=8
hbarom=10

time ./hf << end
&hf_data hbaromega=$hbarom, num_part = $num_part, Nmax = $Nmax, is_RPA = .true., RPA_v_res_scale = 1.0, iter_max = 40 /
end

cp RPA_energies.dat RPA.energies.n$num_part.Nmax$Nmax.hbarom$hbarom.dat
cp RPA_matrix.dat RPA.matrix.n$num_part.Nmax$Nmax.hbarom$hbarom.dat
cp RPA_XY.dat RPA.XY.n$num_part.Nmax$Nmax.hbarom$hbarom.dat
cp RPA_ph_states.dat RPA.phstates.n$num_part.Nmax$Nmax.hbarom$hbarom.dat
cp energies_rpa.agr energies.rpa.n$num_part.Nmax$Nmax.hbarom$hbarom.agr
