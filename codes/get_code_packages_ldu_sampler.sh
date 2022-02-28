#!/usr/bin/env bash

# download the code package

# # download MUSIC
# rm -fr MUSIC_code
# git clone --depth=1 https://github.com/LipeiDu/MUSIC.git -b du_viscous MUSIC_code
# rm -fr MUSIC_code/.git

# # download iSS particle sampler
# rm -fr iSS_code
# git clone --depth=1 https://github.com/chunshen1987/iSS -b momentumSampler iSS_code
# rm -fr iSS_code/.git

# # download iS3D particle sampler
# rm -fr iS3D_code
# git clone --depth=1 https://github.com/LipeiDu/iS3D2.git -b ldu_dev iS3D_code
# rm -fr iS3D_code/.git

# download UrQMD afterburner
rm -fr urqmd_code
git clone --depth=1 https://Chunshen1987@bitbucket.org/Chunshen1987/urqmd_afterburner.git urqmd_code
rm -fr urqmd_code/.git

# download hadronic afterner
# rm -fr hadronic_afterburner_toolkit_code
# git clone --depth=1 https://github.com/chunshen1987/hadronic_afterburner_toolkit hadronic_afterburner_toolkit_code
# rm -fr hadronic_afterburner_toolkit_code/.git

