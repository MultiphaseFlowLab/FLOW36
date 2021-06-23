# clean output, nohup
rm 2D_divergence/nohup.out
rm 2D_divergence/output/*
rm angle_interface_vorticity/nohup.out
rm angle_interface_vorticity/output/*
rm coalescence_breakup/nohup.out
rm coalescence_breakup/output/*
rm curvature_pdf/nohup.out
rm curvature_pdf/output/*
rm curvature_time/nohup.out
rm curvature_time/output/*
rm curvature_time_psi/nohup.out
rm curvature_time_psi/output/*
rm curvature_time_kin/nohup.out
rm curvature_time_kin/output/*
rm drop_count/nohup.out
rm drop_count/output/*
rm eigenvalues_vectors_2D_xyz/nohup.out
rm eigenvalues_vectors_2D_xyz/output/*
rm eigenvalues_vectors_2D_yz/nohup.out
rm eigenvalues_vectors_2D_yz/output/*
rm energy_budget_MKE/nohup.out
rm energy_budget_MKE/output/*
rm generate_pdf/nohup.out
rm generate_pdf/output/*
rm get_interface/nohup.out
rm get_interface/output/*
rm get_pressure_MP/nohup.out
rm get_pressure_MP/output/*
rm invariants/nohup.out
rm invariants/output/*
rm join_2D_div/nohup.out
rm join_2D_div/output/*
rm marangoni/nohup.out
rm marangoni/output/*
rm mass_center/nohup.out
rm mass_center/output/*
rm psi_V/nohup.out
rm psi_V/output/*
rm paraview_finegrid/nohup.out
rm paraview_finegrid/output/*
rm paraview_output_fg/nohup.out
rm paraview_output_fg/output/*
rm paraview_vorticity/nohup.out
rm paraview_vorticity/output/*
rm split_data/nohup.out
rm split_data/output/*
rm stress/nohup.out
rm stress/output/*
rm topology_parameter/nohup.out
rm topology_parameter/output/*
rm wall_shear/nohup.out
rm wall_shear/output/*
rm z_distribution/nohup.out
rm z_distribution/output/*

echo ' '
echo 'sync home s16'
rsync -avz ./* ~/postprocessing_tools/
echo ' '
echo 'sync Marangoni'
rsync -avz ./* marangoni@marangoni.fluid.tuwien.ac.at://home/marangoni/Documents/code/postprocessing_tools/

