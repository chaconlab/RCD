################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../chain_manip.cpp \
../coutsias.cpp \
../ga.cpp \
../matrix.cpp \
../own_io.cpp \
../own_math.cpp \
../rcd_stuff.cpp \
../stat.cpp \
../statistics.cpp 

OBJS += \
./chain_manip.o \
./coutsias.o \
./ga.o \
./matrix.o \
./own_io.o \
./own_math.o \
./rcd_stuff.o \
./stat.o \
./statistics.o 

CPP_DEPS += \
./chain_manip.d \
./coutsias.d \
./ga.d \
./matrix.d \
./own_io.d \
./own_math.d \
./rcd_stuff.d \
./stat.d \
./statistics.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel C++ Compiler'
	icpx -g -O3 -inline-level=2 -I../include -I../../../src -DVOLUME_INCLUDED -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


