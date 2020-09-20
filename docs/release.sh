rsync -av --delete ./build/* leandro@leandro:./public_html/m3g/PDBTools/
ssh leandro lftp -f /home/leandro/programs/scripts/PDBTools_update.lftp

