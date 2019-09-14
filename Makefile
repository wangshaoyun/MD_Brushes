.PHONY: main
main: compt_acc.o glob_var.o main.o \
	in_output.o init_update.o
	
	gfortran -g -Wall -o main compt_acc.o glob_var.o \
	main.o in_output.o init_update.o  -lfftw3
	
# .PHONY: main.o
main.o: initialize_update.mod global_variables.mod \
	input_output.mod compute_acceleration.mod main.f90
	
	gfortran -g -c main.f90
	
# .PHONY: init_update.o initialize_update.mod
init_update.o initialize_update.mod: compute_acceleration.mod global_variables.mod \
	input_output.mod init_update.f90
	
	gfortran -g -c init_update.f90
	
# .PHONY: in_output.o input_output.mod
in_output.o input_output.mod: global_variables.mod compute_acceleration.mod in_output.f90
	
	gfortran -g -c in_output.f90

# .PHONY: compt_acc.o compute_acceleration.mod
compt_acc.o compute_acceleration.mod: global_variables.mod compt_acc.f90
	
	gfortran -g -c compt_acc.f90
	
# .PHONY: glob_var.o global_variables.mod
glob_var.o global_variables.mod: glob_var.f90
	
	gfortran -g -c glob_var.f90
	
clean:
	
	rm compt_acc.o glob_var.o main.o \
		 in_output.o init_update.o \
		 compute_acceleration.mod global_variables.mod \
		 input_output.mod initialize_update.mod \
		 main