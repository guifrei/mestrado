################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F90_SRCS += \
../edge_class_module.f90 \
../element_class_module.f90 \
../linked_list_class_module.f90 \
../main.f90 \
../node_class_module.f90 

OBJS += \
./edge_class_module.o \
./element_class_module.o \
./linked_list_class_module.o \
./main.o \
./node_class_module.o 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -funderscoring -O0 -g -Wall -c -fmessage-length=0 -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

edge_class_module.o: ../edge_class_module.f90 linked_list_class_module.o

element_class_module.o: ../element_class_module.f90 linked_list_class_module.o

linked_list_class_module.o: ../linked_list_class_module.f90

main.o: ../main.f90 edge_class_module.o element_class_module.o linked_list_class_module.o node_class_module.o

node_class_module.o: ../node_class_module.f90 linked_list_class_module.o


