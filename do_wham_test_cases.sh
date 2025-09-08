python do_wham.py windows_reus1_T310_N17_alpha.dat do_wham_py_test -print_progress -plot

python do_wham.py windows_reus1_T310_N17_alpha.dat do_wham_py_test -print_progress -plot -n_part 2
python do_wham.py windows_reus1_T310_N17_alpha.dat do_wham_py_test -print_progress -plot -n_part 3
python do_wham.py windows_reus1_T310_N17_alpha.dat do_wham_py_test -print_progress -plot -n_part 5
python do_wham.py windows_reus1_T310_N17_alpha.dat do_wham_py_test -print_progress -plot -n_part 10

python do_wham.py windows_reus1_T310_N17_alpha.dat do_wham_py_test -t_min 500 -print_progress -plot
python do_wham.py windows_reus1_T310_N17_alpha.dat do_wham_py_test -t_max 500 -print_progress -plot
python do_wham.py windows_reus1_T310_N17_alpha.dat do_wham_py_test -t_min 500 -t_max 400 -print_progress -plot
python do_wham.py windows_reus1_T310_N17_alpha.dat do_wham_py_test -t_min 500 -t_max 500.1 -print_progress -plot
python do_wham.py windows_reus1_T310_N17_alpha.dat do_wham_py_test -t_min 500 -t_max 500.01 -print_progress -plot
python do_wham.py windows_reus1_T310_N17_alpha.dat do_wham_py_test -t_min 500 -t_max 500 -print_progress -plot

python do_wham.py windows_reus1_T310_N17_alpha.dat do_wham_py_test -xi_min 5 -print_progress -plot
python do_wham.py windows_reus1_T310_N17_alpha.dat do_wham_py_test -xi_max 1 -print_progress -plot
python do_wham.py windows_reus1_T310_N17_alpha.dat do_wham_py_test -xi_min 2 -xi_max 3.5 -print_progress -plot
python do_wham.py windows_reus1_T310_N17_alpha.dat do_wham_py_test -xi_min 2 -xi_max 1.9 -print_progress -plot
python do_wham.py windows_reus1_T310_N17_alpha.dat do_wham_py_test -xi_min 2 -xi_max 2.0001 -print_progress -plot
python do_wham.py windows_reus1_T310_N17_alpha.dat do_wham_py_test -xi_min 2 -xi_max 2.00001 -print_progress -plot
python do_wham.py windows_reus1_T310_N17_alpha.dat do_wham_py_test -xi_min 2 -xi_max 2 -print_progress -plot

python do_wham.py windows_reus1_T310_N17_alpha.dat do_wham_py_test -t_min 1000 -xi_max 3.5 -print_progress -plot

# files not located within wham_v7 directory
# python do_wham.py umbrella_test.win.txt do_wham_py_test_flatdata -print_progress -plot -kT 0.5961573
# python do_wham.py umbrella_test.win.txt do_wham_py_test_flatdata -print_progress -plot -kT 0.8
# python do_wham.py umbrella_test.win.txt do_wham_py_test_flatdata -print_progress -plot -kT 0.4

# files not located within wham_v7 directory
# python do_wham.py windows_fulldata_reus1_T310_N17_alpha.dat do_wham_py_test_fulldata -print_progress -plot
