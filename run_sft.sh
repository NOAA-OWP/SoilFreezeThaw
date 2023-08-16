#!/bin/bash
option=$1

if [ ! $# == 1 ]; then
    echo "Usage: $0 OPTION={STANDALONE,PFRAMEWORK}"
    echo "One of these options must be specified to run SFT model"
    exit
fi

if [ $option == "STANDALONE" ] || [ $option == "PFRAMEWORK" ]; then
    echo "SFT model running with option $option"
else
    echo "Invalid option! $option"
    exit
fi


args=" "
exe_name=" "
if [ $option == "STANDALONE" ]; then
    args='./configs/laramie_config_standalone.txt'
    exe_name='sft_standalone'
else if [ $option == "PFRAMEWORK" ]; then
	 args="./configs/laramie_config_cfe.txt ./configs/laramie_config_aorc.txt ./configs/laramie_config_pet.txt ./configs/laramie_config_sft.txt ./configs/laramie_config_smp.txt"
	 exe_name='sft_pframework'
     fi
fi
echo "config file: $args"
./build/${exe_name} $args
