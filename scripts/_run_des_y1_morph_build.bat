@echo off
cd /d "c:\Users\lenovo\Desktop\Bhardwajetal_2024_nature_inclination_angle-main"
set PYTHONIOENCODING=utf-8
set PYTHONUNBUFFERED=1
echo ===== DES Y1 morph RANGE-READ sample %DATE% %TIME% =====>> DES_y1_morph_build.log
python -u scripts/build_des_y1_morph_sample.py --target-rows 500000 --block 5000 >> DES_y1_morph_build.log 2>> DES_y1_morph_build.err
echo ===== DONE %DATE% %TIME% (exit %ERRORLEVEL%) =====>> DES_y1_morph_build.log
