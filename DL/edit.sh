#!/bin/bash

QPs=(22)
Dir=Half

for qp in "${QPs[@]}";
do
	cd ~/git-repos/HM16.9/DL/$Dir/$qp
    List=(emb* lins*-weight.csv lins*-bias.csv bns*-weight.csv bns*-bias.csv bns*-running_mean.csv bns*-running_var.csv mapper*)
	
	for file in *;
	do
        truncate -s-2 $file
        sed -i -e 's/^/\t\t\t/' $file
        echo ";" >> "$file"
	done
	
	for i in "${!List[@]}";
	do
		mv "${List[$i]}" "$((i+1)).${List[$i]}"
	done
done
