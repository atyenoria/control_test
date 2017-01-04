################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../estimation.cpp \
../gnuplot.cpp \
../histogram.cpp \
../hrandom++.cpp \
../lqr_test2.cpp 

O_SRCS += \
../estimation.o \
../gnuplot.o \
../hrandom++.o \
../lqr_test2.o 

OBJS += \
./estimation.o \
./gnuplot.o \
./histogram.o \
./hrandom++.o \
./lqr_test2.o 

CPP_DEPS += \
./estimation.d \
./gnuplot.d \
./histogram.d \
./hrandom++.d \
./lqr_test2.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


