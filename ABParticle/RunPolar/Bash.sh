python3 create_matrix.py;
bash run.sh;
savedir="../data/test06/"
src1dir=$savedir"RunPolar"
mkdir $src1dir
cp -r ../src $savedir
cp -r ./ $src1dir
# cp run_cuda.pl $src1dir
# cp run.sh $src1dir
# cp bacteriaPolar.cu $src1dir
# cp create_matrix.py $src1dir
# # cp draw_polar.py $src1dir
# # cp DrawParticle.py $src1dir
# cp DrawTurbulence.py $src1dir
# cp visual_video.py $src1dir
CUDA_VISIBLE_DEVICES=1 perl run_cuda.pl;
# python3 $src1dir"/draw_polar.py"

mv *.dat $savedir
python3 $src1dir"/DrawTurbulence.py"
python3 $src1dir"/visual_video.py"
