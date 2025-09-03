cp -r start.time 3600.00
blockMesh >> log.txt
setExprFields >> log.txt
date >> log.txt
vssFOAM 2>&1 >> /dev/null

date >> log.txt

