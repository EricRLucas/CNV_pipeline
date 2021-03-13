
species_manifest (){ 
	oldIFS=$IFS
	IFS=","
	read names
	while read $names; do
		echo -e ${sample_id}\\t${species_gambcolu_arabiensis/_/}
	done
	IFS=$oldIFS
}

