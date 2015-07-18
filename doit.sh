export NR_TMPDIR=/disk1/tmp
nohup nice -10 time ./nr -v -v /data/blastdb/gbtrans.faa >/disk1/tmp/gbtrans.nr 2>/disk1/tmp/gbtrans.log &

