



function CODING_POTENTIAL_CALCULATION(){

			./CPC2.py -i $1 -o coding_potential_CPC2.result

			python CPPred.py -i $1 \
				-hex ../Hexamer/Human_Hexamer.tsv \
				-r ../Human_Model/Human.range \
				-mol ../Human_Model/Human.model \
				-spe Human -o coding_potential_CPPred.result

}





 ##this function need a $1 possible_lncrna.fa