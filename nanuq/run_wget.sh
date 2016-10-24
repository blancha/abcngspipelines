#PBS -A feb-684-ab
cd $PBS_O_WORKDIR

echo "Please be advised that the download of archived files can take up to one hour before starting."
J_USER="ablanche"
J_PASS="John.Craig.Venter"
echo -n "j_username=$J_USER&j_password=$J_PASS" > .tmp598745661.dat	
chmod "600" .tmp598745661.dat
if [ "$(which wget)" ]
then   
    echo "Downloading with Wget"   
    wget --no-cookies -c -i readSetLinks.txt --post-file .tmp598745661.dat	
    rm .tmp598745661.dat
else 
    echo "Downloading with Curl"
    xargs -n1 curl -C -J -L -d "@.tmp598745661.dat" -O < readSetLinks.txt 
    rm .tmp598745661.dat
fi
