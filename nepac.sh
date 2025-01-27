date=$(date --date '5 hours ago' '+%Y%m%d')
temphour=$(($(date --date '5 hours ago' '+%H') - $(($(date --date '5 hours ago' '+%H')%6))))

printf -v hour '%02d' $temphour

echo $date $hour

grads -blcx "run script.gs $date $hour"
source ~/anaconda3/etc/profile.d/conda.sh
conda activate py3
python nepac.py $date $hour
