#!/bin/bash/
#SET UP THE DIRECTORIES :
#This is where the bundles are :
FOLDERS="rna_structure/test_bundles/*"
#This is where you want them to go after processing :
mkdir /Users/aniamarczynska/DataspellProjects/dsProject/rna_structure/processed_bundles/
OUTPUTDIR="/Users/aniamarczynska/DataspellProjects/dsProject/rna_structure/processed_bundles/"
mkdir "$OUTPUTDIR"unmatched_chains/
for bundlefol in $FOLDERS
do
  echo "${bundlefol}"/*.txt #e.g. 4v89-pdb-bundle/4v89-chain-id-mapping.txt OK
	pdbid=$(basename -s .txt "$bundlefol"| cut -c1-4)
	echo "$pdbid"
	nbund=$(find "$bundlefol" -type f -name "*.pdb" | wc -l)
	echo "PDB files in bundlefolder $nbund"
	if [ $nbund -gt 1 ]
	then
	  #In this instance there's >1 .pdb file
	  #Using the chain mapping file compare the new and original Chain IDs
		cols=$(tail -n +2 "${bundlefol}"/*.txt | grep -v ^"$pdbid" | awk '{$1=$1};NF' | awk '{OFS=","; $1=$1}1')
		while IFS="," read newchain originalchain # Read lines splitting on commas
        do
                #compare each new chain and the original in the file
                if [ "$newchain" = "$originalchain" ]
                then
                  echo "Equal"
                  unmatched=0
                else
                  #As soon as there's a difference, the loop breaks
                  echo "NOT equal."
                  unmatched=1
                  break
                fi
        done <<< "$cols"
		if [ "$unmatched" = 0 ] #If all chain IDs are the same in the new and original the value is 0
		then
		  echo "Matched chains in $pdbid"
		  #concat the .pdb files
		  cat "${bundlefol}"/*1.pdb > "${OUTPUTDIR}/${pdbid}_bunc.pdb" #first .pdb
		  others=$(find "${bundlefol}" -maxdepth 1 -iname '*.pdb' -not -name '*1.pdb' -exec ls {} +) #for files after ...1.pdb
		  for pdb in $others
		  do
		    echo "$pdb"
		    cat "$pdb" | sed -n '/^ATOM/,$p' >> "${OUTPUTDIR}/${pdbid}_bunc.pdb" #| awk '/^ATOM/ {exit} {print}' #> "${pdbid}_bunc.pdb" #sed -n '1,/^ATOM/!p' > "${pdbid}_bunc.pdb" sed '/^ATOM/q'
		  done
		else
		  echo "unmatched $pdbid"
		  #copy them into directory /unequal_chains
		  cp -r "$bundlefol" "${OUTPUTDIR}unmatched_chains/"
		fi
  else
    #In this instance there's only one .pdb in the bundle folder.
    #Copy .pdb file with new name "[pdbid]...bun.pdb"
    cp "${bundlefol}"/*.pdb "${OUTPUTDIR}/${pdbid}_bun.pdb"
  fi
done
