export NR_TMPDIR=/disk1/tmp
nohup nice -10 time ./nr -n -v -v /disk1/tmp/gbtrans.nr /data/blastdb/pdb260500.faa >/disk1/tmp/gbtrans_pdb.nr 2>/disk1/tmp/gbtrans_pdb.log &

