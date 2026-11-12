import{_ as a,c as n,a5 as i,o as p}from"./chunks/framework.BmjlC9hS.js";const l="/Comrade.jl/v0.11.4/assets/ClosureImaging-42.DD-h7-6H.png",e="/Comrade.jl/v0.11.4/assets/ClosureImaging-44.CGlGA3yY.png",t="/Comrade.jl/v0.11.4/assets/ClosureImaging-52.ChYvhKn9.png",h="/Comrade.jl/v0.11.4/assets/ClosureImaging-54.PibdiavF.png",u=JSON.parse('{"title":"Imaging a Black Hole using only Closure Quantities","description":"","frontmatter":{},"headers":[],"relativePath":"tutorials/intermediate/ClosureImaging.md","filePath":"tutorials/intermediate/ClosureImaging.md","lastUpdated":null}'),r={name:"tutorials/intermediate/ClosureImaging.md"};function k(c,s,o,d,g,E){return p(),n("div",null,s[0]||(s[0]=[i(`<h1 id="Imaging-a-Black-Hole-using-only-Closure-Quantities" tabindex="-1">Imaging a Black Hole using only Closure Quantities <a class="header-anchor" href="#Imaging-a-Black-Hole-using-only-Closure-Quantities" aria-label="Permalink to &quot;Imaging a Black Hole using only Closure Quantities {#Imaging-a-Black-Hole-using-only-Closure-Quantities}&quot;">​</a></h1><p>In this tutorial, we will create a preliminary reconstruction of the 2017 M87 data on April 6 using closure-only imaging. This tutorial is a general introduction to closure-only imaging in Comrade. For an introduction to simultaneous image and instrument modeling, see <a href="/Comrade.jl/v0.11.4/tutorials/intermediate/StokesIImaging#Stokes-I-Simultaneous-Image-and-Instrument-Modeling">Stokes I Simultaneous Image and Instrument Modeling</a></p><h2 id="Introduction-to-Closure-Imaging" tabindex="-1">Introduction to Closure Imaging <a class="header-anchor" href="#Introduction-to-Closure-Imaging" aria-label="Permalink to &quot;Introduction to Closure Imaging {#Introduction-to-Closure-Imaging}&quot;">​</a></h2><p>The EHT is one of the highest-resolution telescope ever created. Its resolution is equivalent to roughly tracking a hockey puck on the moon when viewing it from the earth. However, the EHT is also a unique interferometer. First, EHT data is incredibly sparse, the array is formed from only eight geographic locations around the planet. Second, the obseving frequency is much higher than traditional VLBI. Lastly, each site in the array is unique. They have different dishes, recievers, feeds, and electronics. Putting this all together implies that many of the common imaging techniques struggle to fit the EHT data and explore the uncertainty in both the image and instrument. One way to deal with some of these uncertainties is to not directly fit the data but instead fit closure quantities, which are independent of many of the instrumental effects that plague the data. The types of closure quantities are briefly described in <a href="/Comrade.jl/v0.11.4/vlbi_imaging_problem#Introduction-to-the-VLBI-Imaging-Problem">Introduction to the VLBI Imaging Problem</a>.</p><p>In this tutorial, we will do closure-only modeling of M87 to produce a posterior of images of M87.</p><p>To get started, we will load Comrade</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Comrade</span></span></code></pre></div><p>Pyehtim loads eht-imaging using PythonCall this is necessary to load uvfits files currently.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Pyehtim</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>    CondaPkg Found dependencies: /home/runner/.julia/packages/PythonCall/Nr75f/CondaPkg.toml</span></span>
<span class="line"><span>    CondaPkg Found dependencies: /home/runner/.julia/packages/Pyehtim/Za8ac/CondaPkg.toml</span></span>
<span class="line"><span>    CondaPkg Resolving changes</span></span>
<span class="line"><span>             + ehtim (pip)</span></span>
<span class="line"><span>             + libstdcxx-ng</span></span>
<span class="line"><span>             + numpy</span></span>
<span class="line"><span>             + numpy (pip)</span></span>
<span class="line"><span>             + openssl</span></span>
<span class="line"><span>             + pandas</span></span>
<span class="line"><span>             + pynfft</span></span>
<span class="line"><span>             + python</span></span>
<span class="line"><span>             + setuptools (pip)</span></span>
<span class="line"><span>             + uv</span></span>
<span class="line"><span>    CondaPkg Creating environment</span></span>
<span class="line"><span>             │ /home/runner/.julia/artifacts/7973f2c7725e2d0eef7a95159454c4145f0945a2/bin/micromamba</span></span>
<span class="line"><span>             │ -r /home/runner/.julia/scratchspaces/0b3b1443-0f03-428d-bdfb-f27f9c1191ea/root</span></span>
<span class="line"><span>             │ create</span></span>
<span class="line"><span>             │ -y</span></span>
<span class="line"><span>             │ -p /home/runner/work/Comrade.jl/Comrade.jl/examples/intermediate/ClosureImaging/.CondaPkg/env</span></span>
<span class="line"><span>             │ --override-channels</span></span>
<span class="line"><span>             │ --no-channel-priority</span></span>
<span class="line"><span>             │ libstdcxx-ng[version=&#39;&gt;=3.4,&lt;13.0&#39;]</span></span>
<span class="line"><span>             │ numpy[version=&#39;&gt;=1.24, &lt;2.0&#39;]</span></span>
<span class="line"><span>             │ openssl[version=&#39;&gt;=3, &lt;3.1&#39;]</span></span>
<span class="line"><span>             │ pandas[version=&#39;&lt;2&#39;]</span></span>
<span class="line"><span>             │ pynfft[version=&#39;&gt;1.0&#39;]</span></span>
<span class="line"><span>             │ python[version=&#39;&gt;=3.8,&lt;4&#39;,channel=&#39;conda-forge&#39;,build=&#39;*cpython*&#39;]</span></span>
<span class="line"><span>             │ python[version=&#39;&gt;=3.6,&lt;=3.10&#39;]</span></span>
<span class="line"><span>             │ uv[version=&#39;&gt;=0.4&#39;]</span></span>
<span class="line"><span>             └ -c conda-forge</span></span>
<span class="line"><span>conda-forge/linux-64                                        Using cache</span></span>
<span class="line"><span>conda-forge/noarch                                          Using cache</span></span>
<span class="line"><span></span></span>
<span class="line"><span>Transaction</span></span>
<span class="line"><span></span></span>
<span class="line"><span>  Prefix: /home/runner/work/Comrade.jl/Comrade.jl/examples/intermediate/ClosureImaging/.CondaPkg/env</span></span>
<span class="line"><span></span></span>
<span class="line"><span>  Updating specs:</span></span>
<span class="line"><span></span></span>
<span class="line"><span>   - libstdcxx-ng[version=&#39;&gt;=3.4,&lt;13.0&#39;]</span></span>
<span class="line"><span>   - numpy[version=&#39;&gt;=1.24, &lt;2.0&#39;]</span></span>
<span class="line"><span>   - openssl[version=&#39;&gt;=3, &lt;3.1&#39;]</span></span>
<span class="line"><span>   - pandas[version=&#39;&lt;2&#39;]</span></span>
<span class="line"><span>   - pynfft[version=&#39;&gt;1.0&#39;]</span></span>
<span class="line"><span>   - conda-forge::python[version=&#39;&gt;=3.8,&lt;4&#39;,build=*cpython*]</span></span>
<span class="line"><span>   - python[version=&#39;&gt;=3.6,&lt;=3.10&#39;]</span></span>
<span class="line"><span>   - uv[version=&#39;&gt;=0.4&#39;]</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>  Package                 Version  Build                Channel           Size</span></span>
<span class="line"><span>────────────────────────────────────────────────────────────────────────────────</span></span>
<span class="line"><span>  Install:</span></span>
<span class="line"><span>────────────────────────────────────────────────────────────────────────────────</span></span>
<span class="line"><span></span></span>
<span class="line"><span>  + libstdcxx-ng           12.3.0  hc0a3c3a_7           conda-forge     Cached</span></span>
<span class="line"><span>  + _libgcc_mutex             0.1  conda_forge          conda-forge     Cached</span></span>
<span class="line"><span>  + python_abi               3.10  5_cp310              conda-forge     Cached</span></span>
<span class="line"><span>  + ca-certificates    2024.12.14  hbcca054_0           conda-forge     Cached</span></span>
<span class="line"><span>  + ld_impl_linux-64         2.43  h712a8e2_2           conda-forge     Cached</span></span>
<span class="line"><span>  + libgomp                14.2.0  h77fa898_1           conda-forge     Cached</span></span>
<span class="line"><span>  + _openmp_mutex             4.5  2_gnu                conda-forge     Cached</span></span>
<span class="line"><span>  + libgcc                 14.2.0  h77fa898_1           conda-forge     Cached</span></span>
<span class="line"><span>  + liblzma                 5.6.3  hb9d3cd8_1           conda-forge     Cached</span></span>
<span class="line"><span>  + libzlib                 1.3.1  hb9d3cd8_2           conda-forge     Cached</span></span>
<span class="line"><span>  + libgfortran5           14.2.0  hd5240d6_1           conda-forge     Cached</span></span>
<span class="line"><span>  + libstdcxx              14.2.0  hc0a3c3a_1           conda-forge     Cached</span></span>
<span class="line"><span>  + libgcc-ng              14.2.0  h69a702a_1           conda-forge     Cached</span></span>
<span class="line"><span>  + xz-tools                5.6.3  hb9d3cd8_1           conda-forge     Cached</span></span>
<span class="line"><span>  + xz-gpl-tools            5.6.3  hbcc6ac9_1           conda-forge     Cached</span></span>
<span class="line"><span>  + liblzma-devel           5.6.3  hb9d3cd8_1           conda-forge     Cached</span></span>
<span class="line"><span>  + libsqlite              3.47.2  hee588c1_0           conda-forge     Cached</span></span>
<span class="line"><span>  + libgfortran            14.2.0  h69a702a_1           conda-forge     Cached</span></span>
<span class="line"><span>  + uv                     0.5.11  h0f3a69f_0           conda-forge     Cached</span></span>
<span class="line"><span>  + tk                     8.6.13  noxft_h4845f30_101   conda-forge     Cached</span></span>
<span class="line"><span>  + ncurses                   6.5  he02047a_1           conda-forge     Cached</span></span>
<span class="line"><span>  + libuuid                2.38.1  h0b41bf4_0           conda-forge     Cached</span></span>
<span class="line"><span>  + libnsl                  2.0.1  hd590300_0           conda-forge     Cached</span></span>
<span class="line"><span>  + libffi                  3.4.2  h7f98852_5           conda-forge     Cached</span></span>
<span class="line"><span>  + bzip2                   1.0.8  h4bc722e_7           conda-forge     Cached</span></span>
<span class="line"><span>  + openssl                3.0.14  h4ab18f5_0           conda-forge     Cached</span></span>
<span class="line"><span>  + xz                      5.6.3  hbcc6ac9_1           conda-forge     Cached</span></span>
<span class="line"><span>  + libgfortran-ng         14.2.0  h69a702a_1           conda-forge     Cached</span></span>
<span class="line"><span>  + libopenblas            0.3.28  pthreads_h94d23a6_1  conda-forge     Cached</span></span>
<span class="line"><span>  + readline                  8.2  h8228510_1           conda-forge     Cached</span></span>
<span class="line"><span>  + fftw                   3.3.10  nompi_hf1063bd_110   conda-forge     Cached</span></span>
<span class="line"><span>  + libblas                 3.9.0  26_linux64_openblas  conda-forge     Cached</span></span>
<span class="line"><span>  + sqlite                 3.47.2  h9eae976_0           conda-forge     Cached</span></span>
<span class="line"><span>  + nfft                    3.2.4  hf8c457e_1000        conda-forge     Cached</span></span>
<span class="line"><span>  + libcblas                3.9.0  26_linux64_openblas  conda-forge     Cached</span></span>
<span class="line"><span>  + liblapack               3.9.0  26_linux64_openblas  conda-forge     Cached</span></span>
<span class="line"><span>  + tzdata                  2024b  hc8b5060_0           conda-forge     Cached</span></span>
<span class="line"><span>  + python                 3.10.0  h543edf9_3_cpython   conda-forge     Cached</span></span>
<span class="line"><span>  + wheel                  0.45.1  pyhd8ed1ab_1         conda-forge     Cached</span></span>
<span class="line"><span>  + setuptools             75.6.0  pyhff2d567_1         conda-forge     Cached</span></span>
<span class="line"><span>  + pip                    24.3.1  pyh8b19718_2         conda-forge     Cached</span></span>
<span class="line"><span>  + six                    1.17.0  pyhd8ed1ab_0         conda-forge     Cached</span></span>
<span class="line"><span>  + pytz                   2024.2  pyhd8ed1ab_1         conda-forge     Cached</span></span>
<span class="line"><span>  + python-dateutil   2.9.0.post0  pyhff2d567_1         conda-forge     Cached</span></span>
<span class="line"><span>  + numpy                  1.26.4  py310hb13e2d6_0      conda-forge     Cached</span></span>
<span class="line"><span>  + pynfft                  1.3.2  py310hde88566_1006   conda-forge     Cached</span></span>
<span class="line"><span>  + pandas                  1.5.3  py310h9b08913_1      conda-forge     Cached</span></span>
<span class="line"><span></span></span>
<span class="line"><span>  Summary:</span></span>
<span class="line"><span></span></span>
<span class="line"><span>  Install: 47 packages</span></span>
<span class="line"><span></span></span>
<span class="line"><span>  Total download: 0 B</span></span>
<span class="line"><span></span></span>
<span class="line"><span>────────────────────────────────────────────────────────────────────────────────</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>Transaction starting</span></span>
<span class="line"><span>Linking libstdcxx-ng-12.3.0-hc0a3c3a_7</span></span>
<span class="line"><span>Linking _libgcc_mutex-0.1-conda_forge</span></span>
<span class="line"><span>Linking python_abi-3.10-5_cp310</span></span>
<span class="line"><span>Linking ca-certificates-2024.12.14-hbcca054_0</span></span>
<span class="line"><span>Linking ld_impl_linux-64-2.43-h712a8e2_2</span></span>
<span class="line"><span>Linking libgomp-14.2.0-h77fa898_1</span></span>
<span class="line"><span>Linking _openmp_mutex-4.5-2_gnu</span></span>
<span class="line"><span>Linking libgcc-14.2.0-h77fa898_1</span></span>
<span class="line"><span>Linking liblzma-5.6.3-hb9d3cd8_1</span></span>
<span class="line"><span>Linking libzlib-1.3.1-hb9d3cd8_2</span></span>
<span class="line"><span>Linking libgfortran5-14.2.0-hd5240d6_1</span></span>
<span class="line"><span>Linking libstdcxx-14.2.0-hc0a3c3a_1</span></span>
<span class="line"><span>warning  libmamba [libstdcxx-14.2.0-hc0a3c3a_1] The following files were already present in the environment:</span></span>
<span class="line"><span>    - lib/libstdc++.so</span></span>
<span class="line"><span>    - lib/libstdc++.so.6</span></span>
<span class="line"><span>    - share/licenses/libstdc++/RUNTIME.LIBRARY.EXCEPTION</span></span>
<span class="line"><span>Linking libgcc-ng-14.2.0-h69a702a_1</span></span>
<span class="line"><span>Linking xz-tools-5.6.3-hb9d3cd8_1</span></span>
<span class="line"><span>Linking xz-gpl-tools-5.6.3-hbcc6ac9_1</span></span>
<span class="line"><span>Linking liblzma-devel-5.6.3-hb9d3cd8_1</span></span>
<span class="line"><span>Linking libsqlite-3.47.2-hee588c1_0</span></span>
<span class="line"><span>Linking libgfortran-14.2.0-h69a702a_1</span></span>
<span class="line"><span>Linking uv-0.5.11-h0f3a69f_0</span></span>
<span class="line"><span>Linking tk-8.6.13-noxft_h4845f30_101</span></span>
<span class="line"><span>Linking ncurses-6.5-he02047a_1</span></span>
<span class="line"><span>Linking libuuid-2.38.1-h0b41bf4_0</span></span>
<span class="line"><span>Linking libnsl-2.0.1-hd590300_0</span></span>
<span class="line"><span>Linking libffi-3.4.2-h7f98852_5</span></span>
<span class="line"><span>Linking bzip2-1.0.8-h4bc722e_7</span></span>
<span class="line"><span>Linking openssl-3.0.14-h4ab18f5_0</span></span>
<span class="line"><span>Linking xz-5.6.3-hbcc6ac9_1</span></span>
<span class="line"><span>Linking libgfortran-ng-14.2.0-h69a702a_1</span></span>
<span class="line"><span>Linking libopenblas-0.3.28-pthreads_h94d23a6_1</span></span>
<span class="line"><span>Linking readline-8.2-h8228510_1</span></span>
<span class="line"><span>Linking fftw-3.3.10-nompi_hf1063bd_110</span></span>
<span class="line"><span>Linking libblas-3.9.0-26_linux64_openblas</span></span>
<span class="line"><span>Linking sqlite-3.47.2-h9eae976_0</span></span>
<span class="line"><span>Linking nfft-3.2.4-hf8c457e_1000</span></span>
<span class="line"><span>Linking libcblas-3.9.0-26_linux64_openblas</span></span>
<span class="line"><span>Linking liblapack-3.9.0-26_linux64_openblas</span></span>
<span class="line"><span>Linking tzdata-2024b-hc8b5060_0</span></span>
<span class="line"><span>Linking python-3.10.0-h543edf9_3_cpython</span></span>
<span class="line"><span>Linking wheel-0.45.1-pyhd8ed1ab_1</span></span>
<span class="line"><span>Linking setuptools-75.6.0-pyhff2d567_1</span></span>
<span class="line"><span>Linking pip-24.3.1-pyh8b19718_2</span></span>
<span class="line"><span>Linking six-1.17.0-pyhd8ed1ab_0</span></span>
<span class="line"><span>Linking pytz-2024.2-pyhd8ed1ab_1</span></span>
<span class="line"><span>Linking python-dateutil-2.9.0.post0-pyhff2d567_1</span></span>
<span class="line"><span>Linking numpy-1.26.4-py310hb13e2d6_0</span></span>
<span class="line"><span>Linking pynfft-1.3.2-py310hde88566_1006</span></span>
<span class="line"><span>Linking pandas-1.5.3-py310h9b08913_1</span></span>
<span class="line"><span></span></span>
<span class="line"><span>Transaction finished</span></span>
<span class="line"><span></span></span>
<span class="line"><span>To activate this environment, use:</span></span>
<span class="line"><span></span></span>
<span class="line"><span>    micromamba activate /home/runner/work/Comrade.jl/Comrade.jl/examples/intermediate/ClosureImaging/.CondaPkg/env</span></span>
<span class="line"><span></span></span>
<span class="line"><span>Or to execute a single command in this environment, use:</span></span>
<span class="line"><span></span></span>
<span class="line"><span>    micromamba run -p /home/runner/work/Comrade.jl/Comrade.jl/examples/intermediate/ClosureImaging/.CondaPkg/env mycommand</span></span>
<span class="line"><span></span></span>
<span class="line"><span>    CondaPkg Installing Pip packages</span></span>
<span class="line"><span>             │ /home/runner/work/Comrade.jl/Comrade.jl/examples/intermediate/ClosureImaging/.CondaPkg/env/bin/uv</span></span>
<span class="line"><span>             │ pip</span></span>
<span class="line"><span>             │ install</span></span>
<span class="line"><span>             │ ehtim &gt;=1.2.9, &lt;2.0</span></span>
<span class="line"><span>             │ numpy &gt;=1.24, &lt;2.0</span></span>
<span class="line"><span>             └ setuptools</span></span>
<span class="line"><span>Using Python 3.10.0 environment at: /home/runner/work/Comrade.jl/Comrade.jl/examples/intermediate/ClosureImaging/.CondaPkg/env</span></span>
<span class="line"><span>Resolved 36 packages in 19ms</span></span>
<span class="line"><span>Installed 30 packages in 60ms</span></span>
<span class="line"><span> + astropy==6.1.7</span></span>
<span class="line"><span> + astropy-iers-data==0.2024.12.23.0.33.24</span></span>
<span class="line"><span> + certifi==2024.12.14</span></span>
<span class="line"><span> + charset-normalizer==3.4.0</span></span>
<span class="line"><span> + contourpy==1.3.1</span></span>
<span class="line"><span> + cycler==0.12.1</span></span>
<span class="line"><span> + ehtim==1.2.9</span></span>
<span class="line"><span> + fonttools==4.55.3</span></span>
<span class="line"><span> + future==1.0.0</span></span>
<span class="line"><span> + h5py==3.12.1</span></span>
<span class="line"><span> + hdrhistogram==0.10.3</span></span>
<span class="line"><span> + idna==3.10</span></span>
<span class="line"><span> + jplephem==2.22</span></span>
<span class="line"><span> + kiwisolver==1.4.7</span></span>
<span class="line"><span> + matplotlib==3.10.0</span></span>
<span class="line"><span> + networkx==3.4.2</span></span>
<span class="line"><span> + packaging==24.2</span></span>
<span class="line"><span> + pandas-appender==0.9.8.4</span></span>
<span class="line"><span> + paramsurvey==0.4.20</span></span>
<span class="line"><span> + pbr==6.1.0</span></span>
<span class="line"><span> + pillow==11.0.0</span></span>
<span class="line"><span> + psutil==6.1.1</span></span>
<span class="line"><span> + pyerfa==2.0.1.5</span></span>
<span class="line"><span> + pyparsing==3.2.0</span></span>
<span class="line"><span> + pyyaml==6.0.2</span></span>
<span class="line"><span> + requests==2.32.3</span></span>
<span class="line"><span> + scipy==1.13.1</span></span>
<span class="line"><span> + sgp4==2.23</span></span>
<span class="line"><span> + skyfield==1.49</span></span>
<span class="line"><span> + urllib3==2.3.0</span></span></code></pre></div><p>For reproducibility we use a stable random number genreator</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> StableRNGs</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">rng </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> StableRNG</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">123</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>StableRNGs.LehmerRNG(state=0x000000000000000000000000000000f7)</span></span></code></pre></div><h2 id="Load-the-Data" tabindex="-1">Load the Data <a class="header-anchor" href="#Load-the-Data" aria-label="Permalink to &quot;Load the Data {#Load-the-Data}&quot;">​</a></h2><p>To download the data visit <a href="https://doi.org/10.25739/g85n-f134" target="_blank" rel="noreferrer">https://doi.org/10.25739/g85n-f134</a> To load the eht-imaging obsdata object we do:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">obs </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> ehtim</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">obsdata</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">load_uvfits</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">joinpath</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(__DIR, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;..&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;..&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;Data&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;SR1_M87_2017_096_lo_hops_netcal_StokesI.uvfits&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">))</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>Python: &lt;ehtim.obsdata.Obsdata object at 0x7f1ad6b20190&gt;</span></span></code></pre></div><p>Now we do some minor preprocessing:</p><ul><li><p>Scan average the data since the data have been preprocessed so that the gain phases are coherent.</p></li><li><p>Add 2% systematic noise to deal with calibration issues such as leakage.</p></li></ul><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">obs </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> scan_average</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(obs)</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">add_fractional_noise</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.02</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>Python: &lt;ehtim.obsdata.Obsdata object at 0x7f1ad7985960&gt;</span></span></code></pre></div><p>Now, we extract our closure quantities from the EHT data set. We flag now SNR points since the closure likelihood we use is only applicable to high SNR data.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">dlcamp, dcphase  </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> extract_table</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(obs, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">LogClosureAmplitudes</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(;snrcut</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">3</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">), </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">ClosurePhases</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(;snrcut</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">3</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">))</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>(EHTObservationTable{Comrade.EHTLogClosureAmplitudeDatum{:I}}</span></span>
<span class="line"><span>  source:      M87</span></span>
<span class="line"><span>  mjd:         57849</span></span>
<span class="line"><span>  bandwidth:   1.856e9</span></span>
<span class="line"><span>  sites:       [:AA, :AP, :AZ, :JC, :LM, :PV, :SM]</span></span>
<span class="line"><span>  nsamples:    128, EHTObservationTable{Comrade.EHTClosurePhaseDatum{:I}}</span></span>
<span class="line"><span>  source:      M87</span></span>
<span class="line"><span>  mjd:         57849</span></span>
<span class="line"><span>  bandwidth:   1.856e9</span></span>
<span class="line"><span>  sites:       [:AA, :AP, :AZ, :JC, :LM, :PV, :SM]</span></span>
<span class="line"><span>  nsamples:    152)</span></span></code></pre></div><div class="tip custom-block"><p class="custom-block-title">Note</p><p>Fitting low SNR closure data is complicated and requires a more sophisticated likelihood. If low-SNR data is very important we recommend fitting visibilties with a instrumental model.</p></div><h2 id="Build-the-Model/Posterior" tabindex="-1">Build the Model/Posterior <a class="header-anchor" href="#Build-the-Model/Posterior" aria-label="Permalink to &quot;Build the Model/Posterior {#Build-the-Model/Posterior}&quot;">​</a></h2><p>For our model, we will be using an image model that consists of a raster of point sources, convolved with some pulse or kernel to make a <code>ContinuousImage</code>. To define this model we define the standard two argument function <code>sky</code> that defines the sky model we want to fit. The first argument are the model parameters, and are typically a NamedTuple. The second argument defines the metadata for the model that is typically constant. For our model the constant <code>metdata</code> will just by the mean or prior image.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">function</span><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;"> sky</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(θ, metadata)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    (;fg, c, σimg) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> θ</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    (;mimg) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> metadata</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">    # Apply the GMRF fluctuations to the image</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    rast </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> apply_fluctuations</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">CenteredLR</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(), mimg, σimg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">c</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">params)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    m </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> ContinuousImage</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(((</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">fg))</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">rast, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">BSplinePulse{3}</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">())</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">    # Force the image centroid to be at the origin</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    x0, y0 </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> centroid</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(m)</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">    # Add a large-scale gaussian to deal with the over-resolved mas flux</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    g </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> modify</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Gaussian</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(), </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Stretch</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">μas2rad</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">250.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">), </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">μas2rad</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">250.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)), </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Renormalize</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(fg))</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    return</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> shifted</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(m, </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">x0, </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">y0) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">+</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> g</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>sky (generic function with 1 method)</span></span></code></pre></div><p>Now, let&#39;s set up our image model. The EHT&#39;s nominal resolution is 20-25 μas. Additionally, the EHT is not very sensitive to a larger field of views; typically, 60-80 μas is enough to describe the compact flux of M87. Given this, we only need to use a small number of pixels to describe our image.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">npix </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 32</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">fovxy </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> μas2rad</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">150.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>7.27220521664304e-10</span></span></code></pre></div><p>To define the image model we need to specify both the grid we will be using and the FT algorithm we will use, in this case the NFFT which is the most efficient.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">grid </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> imagepixels</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(fovxy, fovxy, npix, npix)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>RectiGrid(</span></span>
<span class="line"><span>executor: ComradeBase.Serial()</span></span>
<span class="line"><span>Dimensions: </span></span>
<span class="line"><span>(↓ X Sampled{Float64} LinRange{Float64}(-3.5224744018114725e-10, 3.5224744018114725e-10, 32) ForwardOrdered Regular Points,</span></span>
<span class="line"><span>→ Y Sampled{Float64} LinRange{Float64}(-3.5224744018114725e-10, 3.5224744018114725e-10, 32) ForwardOrdered Regular Points)</span></span>
<span class="line"><span>)</span></span></code></pre></div><p>Now we need to specify our image prior. For this work we will use a Gaussian Markov Random field prior</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> VLBIImagePriors, Distributions</span></span></code></pre></div><p>Since we are using a Gaussian Markov random field prior we need to first specify our <code>mean</code> image. For this work we will use a symmetric Gaussian with a FWHM of 50 μas</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">fwhmfac </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 2</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">sqrt</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">log</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">))</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">mpr </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> modify</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Gaussian</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(), </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Stretch</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">μas2rad</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">50.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">./</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">fwhmfac))</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">imgpr </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> intensitymap</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(mpr, grid)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">skymeta </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (;mimg </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> imgpr</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">./</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">flux</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(imgpr));</span></span></code></pre></div><p>Now we can finally form our image prior. For this we use a heirarchical prior where the direct log-ratio image prior is a Gaussian Markov Random Field. The correlation length of the GMRF is a hyperparameter that is fit during imaging. We pass the data to the prior to estimate what the maximumal resolutoin of the array is and prevent the prior from allowing correlation lengths that are much small than the telescope beam size. Note that this GMRF prior has unit variance. For more information on the GMRF prior see the <a href="/Comrade.jl/v0.11.4/api#Comrade.corr_image_prior"><code>corr_image_prior</code></a> doc string.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">cprior </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> corr_image_prior</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(grid, dlcamp)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>HierarchicalPrior(</span></span>
<span class="line"><span>	map: </span></span>
<span class="line"><span>	ConditionalMarkov(</span></span>
<span class="line"><span>Random Field: VLBIImagePriors.GaussMarkovRandomField</span></span>
<span class="line"><span>Graph: MarkovRandomFieldGraph{1}(</span></span>
<span class="line"><span>dims: (32, 32)</span></span>
<span class="line"><span>)</span></span>
<span class="line"><span>)	hyper prior: </span></span>
<span class="line"><span>	Truncated(Distributions.InverseGamma{Float64}(</span></span>
<span class="line"><span>invd: Distributions.Gamma{Float64}(α=1.0, θ=0.0407685911951416)</span></span>
<span class="line"><span>θ: 24.528686684644875</span></span>
<span class="line"><span>)</span></span>
<span class="line"><span>; lower=1.0, upper=64.0)</span></span>
<span class="line"><span></span></span>
<span class="line"><span>)</span></span></code></pre></div><p>Putting everything together the total prior is then our image prior, a prior on the standard deviation of the MRF, and a prior on the fractional flux of the Gaussian component.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">prior </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (c </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> cprior, σimg </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Exponential</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">), fg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Uniform</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">))</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>(c = HierarchicalPrior(</span></span>
<span class="line"><span>	map: </span></span>
<span class="line"><span>	ConditionalMarkov(</span></span>
<span class="line"><span>Random Field: VLBIImagePriors.GaussMarkovRandomField</span></span>
<span class="line"><span>Graph: MarkovRandomFieldGraph{1}(</span></span>
<span class="line"><span>dims: (32, 32)</span></span>
<span class="line"><span>)</span></span>
<span class="line"><span>)	hyper prior: </span></span>
<span class="line"><span>	Truncated(Distributions.InverseGamma{Float64}(</span></span>
<span class="line"><span>invd: Distributions.Gamma{Float64}(α=1.0, θ=0.0407685911951416)</span></span>
<span class="line"><span>θ: 24.528686684644875</span></span>
<span class="line"><span>)</span></span>
<span class="line"><span>; lower=1.0, upper=64.0)</span></span>
<span class="line"><span></span></span>
<span class="line"><span>)</span></span>
<span class="line"><span>, σimg = Distributions.Exponential{Float64}(θ=0.1), fg = Distributions.Uniform{Float64}(a=0.0, b=1.0))</span></span></code></pre></div><p>We can then define our sky model.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">skym </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> SkyModel</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(sky, prior, grid; metadata</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">skymeta)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>SkyModel</span></span>
<span class="line"><span>  with map: sky</span></span>
<span class="line"><span>   on grid: ComradeBase.RectiGrid</span></span></code></pre></div><p>Since we are fitting closures we do not need to include an instrument model, since the closure likelihood is approximately independent of gains in the high SNR limit.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Enzyme</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">post </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> VLBIPosterior</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(skym, dlcamp, dcphase; admode</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">set_runtime_activity</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(Enzyme</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Reverse))</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>VLBIPosterior</span></span>
<span class="line"><span>ObservedSkyModel</span></span>
<span class="line"><span>  with map: sky</span></span>
<span class="line"><span>   on grid: VLBISkyModels.FourierDualDomainIdealInstrumentModelData Products: Comrade.EHTLogClosureAmplitudeDatumComrade.EHTClosurePhaseDatum</span></span></code></pre></div><h2 id="Reconstructing-the-Image" tabindex="-1">Reconstructing the Image <a class="header-anchor" href="#Reconstructing-the-Image" aria-label="Permalink to &quot;Reconstructing the Image {#Reconstructing-the-Image}&quot;">​</a></h2><p>To reconstruct the image we will first use the MAP estimate. This is approach is basically a re-implentation of regularized maximum likelihood (RML) imaging. However, unlike traditional RML imaging we also fit the regularizer hyperparameters, thanks to our interpretation of as our imaging prior as a hierarchical model.</p><p>To optimize our posterior <code>Comrade</code> provides the <code>comrade_opt</code> function. To use this functionality a user first needs to import <code>Optimization.jl</code> and the optimizer of choice. In this tutorial we will use Optim.jl&#39;s L-BFGS optimizer, which is defined in the sub-package OptimizationOptimJL. We also need to import Enzyme to allow for automatic differentiation.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Optimization</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> OptimizationOptimJL</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">xopt, sol </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> comrade_opt</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(post, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">LBFGS</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">();</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">                        maxiters</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1000</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, initial_params</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">prior_sample</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(rng, post));</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>┌ Warning: Using fallback BLAS replacements for ([&quot;cblas_zdotc_sub64_&quot;]), performance may be degraded</span></span>
<span class="line"><span>└ @ Enzyme.Compiler ~/.julia/packages/GPUCompiler/GnbhK/src/utils.jl:59</span></span>
<span class="line"><span>┌ Warning: TODO forward zero-set of arraycopy used memset rather than runtime type </span></span>
<span class="line"><span>│ Caused by:</span></span>
<span class="line"><span>│ Stacktrace:</span></span>
<span class="line"><span>│  [1] copy</span></span>
<span class="line"><span>│    @ ./array.jl:411</span></span>
<span class="line"><span>│  [2] map</span></span>
<span class="line"><span>│    @ ./tuple.jl:293</span></span>
<span class="line"><span>│  [3] map</span></span>
<span class="line"><span>│    @ ./tuple.jl:294</span></span>
<span class="line"><span>│  [4] map</span></span>
<span class="line"><span>│    @ ./namedtuple.jl:265</span></span>
<span class="line"><span>│  [5] copy</span></span>
<span class="line"><span>│    @ ~/.julia/packages/StructArrays/KLKnN/src/structarray.jl:467</span></span>
<span class="line"><span>└ @ Enzyme.Compiler ~/.julia/packages/GPUCompiler/GnbhK/src/utils.jl:59</span></span>
<span class="line"><span>┌ Warning: TODO forward zero-set of arraycopy used memset rather than runtime type </span></span>
<span class="line"><span>│ Caused by:</span></span>
<span class="line"><span>│ Stacktrace:</span></span>
<span class="line"><span>│  [1] copy</span></span>
<span class="line"><span>│    @ ./array.jl:411</span></span>
<span class="line"><span>│  [2] map</span></span>
<span class="line"><span>│    @ ./tuple.jl:293</span></span>
<span class="line"><span>│  [3] map</span></span>
<span class="line"><span>│    @ ./tuple.jl:294</span></span>
<span class="line"><span>│  [4] map</span></span>
<span class="line"><span>│    @ ./namedtuple.jl:265</span></span>
<span class="line"><span>│  [5] copy</span></span>
<span class="line"><span>│    @ ~/.julia/packages/StructArrays/KLKnN/src/structarray.jl:467</span></span>
<span class="line"><span>└ @ Enzyme.Compiler ~/.julia/packages/GPUCompiler/GnbhK/src/utils.jl:59</span></span>
<span class="line"><span>┌ Warning: TODO forward zero-set of arraycopy used memset rather than runtime type </span></span>
<span class="line"><span>│ Caused by:</span></span>
<span class="line"><span>│ Stacktrace:</span></span>
<span class="line"><span>│  [1] copy</span></span>
<span class="line"><span>│    @ ./array.jl:411</span></span>
<span class="line"><span>│  [2] map</span></span>
<span class="line"><span>│    @ ./tuple.jl:293</span></span>
<span class="line"><span>│  [3] map</span></span>
<span class="line"><span>│    @ ./tuple.jl:294</span></span>
<span class="line"><span>│  [4] map</span></span>
<span class="line"><span>│    @ ./namedtuple.jl:265</span></span>
<span class="line"><span>│  [5] copy</span></span>
<span class="line"><span>│    @ ~/.julia/packages/StructArrays/KLKnN/src/structarray.jl:467</span></span>
<span class="line"><span>└ @ Enzyme.Compiler ~/.julia/packages/GPUCompiler/GnbhK/src/utils.jl:59</span></span>
<span class="line"><span>┌ Warning: TODO forward zero-set of arraycopy used memset rather than runtime type </span></span>
<span class="line"><span>│ Caused by:</span></span>
<span class="line"><span>│ Stacktrace:</span></span>
<span class="line"><span>│  [1] copy</span></span>
<span class="line"><span>│    @ ./array.jl:411</span></span>
<span class="line"><span>│  [2] map</span></span>
<span class="line"><span>│    @ ./tuple.jl:294</span></span>
<span class="line"><span>│  [3] map</span></span>
<span class="line"><span>│    @ ./namedtuple.jl:265</span></span>
<span class="line"><span>│  [4] copy</span></span>
<span class="line"><span>│    @ ~/.julia/packages/StructArrays/KLKnN/src/structarray.jl:467</span></span>
<span class="line"><span>└ @ Enzyme.Compiler ~/.julia/packages/GPUCompiler/GnbhK/src/utils.jl:59</span></span>
<span class="line"><span>┌ Warning: TODO forward zero-set of arraycopy used memset rather than runtime type </span></span>
<span class="line"><span>│ Caused by:</span></span>
<span class="line"><span>│ Stacktrace:</span></span>
<span class="line"><span>│  [1] copy</span></span>
<span class="line"><span>│    @ ./array.jl:411</span></span>
<span class="line"><span>│  [2] map</span></span>
<span class="line"><span>│    @ ./tuple.jl:293</span></span>
<span class="line"><span>│  [3] map</span></span>
<span class="line"><span>│    @ ./tuple.jl:294</span></span>
<span class="line"><span>│  [4] map</span></span>
<span class="line"><span>│    @ ./namedtuple.jl:265</span></span>
<span class="line"><span>│  [5] copy</span></span>
<span class="line"><span>│    @ ~/.julia/packages/StructArrays/KLKnN/src/structarray.jl:467</span></span>
<span class="line"><span>└ @ Enzyme.Compiler ~/.julia/packages/GPUCompiler/GnbhK/src/utils.jl:59</span></span>
<span class="line"><span>┌ Warning: TODO forward zero-set of arraycopy used memset rather than runtime type </span></span>
<span class="line"><span>│ Caused by:</span></span>
<span class="line"><span>│ Stacktrace:</span></span>
<span class="line"><span>│  [1] copy</span></span>
<span class="line"><span>│    @ ./array.jl:411</span></span>
<span class="line"><span>│  [2] map</span></span>
<span class="line"><span>│    @ ./tuple.jl:293</span></span>
<span class="line"><span>│  [3] map</span></span>
<span class="line"><span>│    @ ./tuple.jl:294</span></span>
<span class="line"><span>│  [4] map</span></span>
<span class="line"><span>│    @ ./namedtuple.jl:265</span></span>
<span class="line"><span>│  [5] copy</span></span>
<span class="line"><span>│    @ ~/.julia/packages/StructArrays/KLKnN/src/structarray.jl:467</span></span>
<span class="line"><span>└ @ Enzyme.Compiler ~/.julia/packages/GPUCompiler/GnbhK/src/utils.jl:59</span></span>
<span class="line"><span>┌ Warning: TODO forward zero-set of arraycopy used memset rather than runtime type </span></span>
<span class="line"><span>│ Caused by:</span></span>
<span class="line"><span>│ Stacktrace:</span></span>
<span class="line"><span>│  [1] copy</span></span>
<span class="line"><span>│    @ ./array.jl:411</span></span>
<span class="line"><span>│  [2] map</span></span>
<span class="line"><span>│    @ ./tuple.jl:293</span></span>
<span class="line"><span>│  [3] map</span></span>
<span class="line"><span>│    @ ./tuple.jl:294</span></span>
<span class="line"><span>│  [4] map</span></span>
<span class="line"><span>│    @ ./namedtuple.jl:265</span></span>
<span class="line"><span>│  [5] copy</span></span>
<span class="line"><span>│    @ ~/.julia/packages/StructArrays/KLKnN/src/structarray.jl:467</span></span>
<span class="line"><span>└ @ Enzyme.Compiler ~/.julia/packages/GPUCompiler/GnbhK/src/utils.jl:59</span></span>
<span class="line"><span>┌ Warning: TODO forward zero-set of arraycopy used memset rather than runtime type </span></span>
<span class="line"><span>│ Caused by:</span></span>
<span class="line"><span>│ Stacktrace:</span></span>
<span class="line"><span>│  [1] copy</span></span>
<span class="line"><span>│    @ ./array.jl:411</span></span>
<span class="line"><span>│  [2] map</span></span>
<span class="line"><span>│    @ ./tuple.jl:294</span></span>
<span class="line"><span>│  [3] map</span></span>
<span class="line"><span>│    @ ./namedtuple.jl:265</span></span>
<span class="line"><span>│  [4] copy</span></span>
<span class="line"><span>│    @ ~/.julia/packages/StructArrays/KLKnN/src/structarray.jl:467</span></span>
<span class="line"><span>└ @ Enzyme.Compiler ~/.julia/packages/GPUCompiler/GnbhK/src/utils.jl:59</span></span></code></pre></div><p>First we will evaluate our fit by plotting the residuals</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Plots</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">p </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> residual</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(post, xopt)</span></span></code></pre></div><p><img src="`+l+`" alt=""></p><p>Now let&#39;s plot the MAP estimate.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">import</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> CairoMakie </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">as</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> CM</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">g </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> imagepixels</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">μas2rad</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">150.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">), </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">μas2rad</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">150.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">), </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">100</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">100</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">img </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> intensitymap</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">skymodel</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(post, xopt), g)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">fig </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> imageviz</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(img, size</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">600</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">500</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">));</span></span></code></pre></div><p><img src="`+e+`" alt=""></p><p>That doesn&#39;t look great. This is pretty common for the sparse EHT data. In this case the MAP often drastically overfits the data, producing a image filled with artifacts. In addition, we note that the MAP itself is not invariant to the model parameterization. Namely, if we changed our prior to use a fully centered parameterization we would get a very different image. Fortunately, these issues go away when we sample from the posterior, and construct expectations of the posterior, like the mean image.</p><p>To sample from the posterior we will use HMC and more specifically the NUTS algorithm. For information about NUTS see Michael Betancourt&#39;s <a href="https://arxiv.org/abs/1701.02434" target="_blank" rel="noreferrer">notes</a>.</p><div class="tip custom-block"><p class="custom-block-title">Note</p><p>For our <code>metric</code> we use a diagonal matrix due to easier tuning.</p></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> AdvancedHMC</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">chain </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> sample</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(rng, post, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">NUTS</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.8</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">), </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">700</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">; n_adapts</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">500</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, progress</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">false</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, initial_params</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">xopt);</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>┌ Warning: Using fallback BLAS replacements for ([&quot;cblas_zdotc_sub64_&quot;]), performance may be degraded</span></span>
<span class="line"><span>└ @ Enzyme.Compiler ~/.julia/packages/GPUCompiler/GnbhK/src/utils.jl:59</span></span>
<span class="line"><span>┌ Warning: TODO forward zero-set of arraycopy used memset rather than runtime type </span></span>
<span class="line"><span>│ Caused by:</span></span>
<span class="line"><span>│ Stacktrace:</span></span>
<span class="line"><span>│  [1] copy</span></span>
<span class="line"><span>│    @ ./array.jl:411</span></span>
<span class="line"><span>│  [2] map</span></span>
<span class="line"><span>│    @ ./tuple.jl:293</span></span>
<span class="line"><span>│  [3] map</span></span>
<span class="line"><span>│    @ ./tuple.jl:294</span></span>
<span class="line"><span>│  [4] map</span></span>
<span class="line"><span>│    @ ./namedtuple.jl:265</span></span>
<span class="line"><span>│  [5] copy</span></span>
<span class="line"><span>│    @ ~/.julia/packages/StructArrays/KLKnN/src/structarray.jl:467</span></span>
<span class="line"><span>└ @ Enzyme.Compiler ~/.julia/packages/GPUCompiler/GnbhK/src/utils.jl:59</span></span>
<span class="line"><span>┌ Warning: TODO forward zero-set of arraycopy used memset rather than runtime type </span></span>
<span class="line"><span>│ Caused by:</span></span>
<span class="line"><span>│ Stacktrace:</span></span>
<span class="line"><span>│  [1] copy</span></span>
<span class="line"><span>│    @ ./array.jl:411</span></span>
<span class="line"><span>│  [2] map</span></span>
<span class="line"><span>│    @ ./tuple.jl:293</span></span>
<span class="line"><span>│  [3] map</span></span>
<span class="line"><span>│    @ ./tuple.jl:294</span></span>
<span class="line"><span>│  [4] map</span></span>
<span class="line"><span>│    @ ./namedtuple.jl:265</span></span>
<span class="line"><span>│  [5] copy</span></span>
<span class="line"><span>│    @ ~/.julia/packages/StructArrays/KLKnN/src/structarray.jl:467</span></span>
<span class="line"><span>└ @ Enzyme.Compiler ~/.julia/packages/GPUCompiler/GnbhK/src/utils.jl:59</span></span>
<span class="line"><span>┌ Warning: TODO forward zero-set of arraycopy used memset rather than runtime type </span></span>
<span class="line"><span>│ Caused by:</span></span>
<span class="line"><span>│ Stacktrace:</span></span>
<span class="line"><span>│  [1] copy</span></span>
<span class="line"><span>│    @ ./array.jl:411</span></span>
<span class="line"><span>│  [2] map</span></span>
<span class="line"><span>│    @ ./tuple.jl:293</span></span>
<span class="line"><span>│  [3] map</span></span>
<span class="line"><span>│    @ ./tuple.jl:294</span></span>
<span class="line"><span>│  [4] map</span></span>
<span class="line"><span>│    @ ./namedtuple.jl:265</span></span>
<span class="line"><span>│  [5] copy</span></span>
<span class="line"><span>│    @ ~/.julia/packages/StructArrays/KLKnN/src/structarray.jl:467</span></span>
<span class="line"><span>└ @ Enzyme.Compiler ~/.julia/packages/GPUCompiler/GnbhK/src/utils.jl:59</span></span>
<span class="line"><span>┌ Warning: TODO forward zero-set of arraycopy used memset rather than runtime type </span></span>
<span class="line"><span>│ Caused by:</span></span>
<span class="line"><span>│ Stacktrace:</span></span>
<span class="line"><span>│  [1] copy</span></span>
<span class="line"><span>│    @ ./array.jl:411</span></span>
<span class="line"><span>│  [2] map</span></span>
<span class="line"><span>│    @ ./tuple.jl:294</span></span>
<span class="line"><span>│  [3] map</span></span>
<span class="line"><span>│    @ ./namedtuple.jl:265</span></span>
<span class="line"><span>│  [4] copy</span></span>
<span class="line"><span>│    @ ~/.julia/packages/StructArrays/KLKnN/src/structarray.jl:467</span></span>
<span class="line"><span>└ @ Enzyme.Compiler ~/.julia/packages/GPUCompiler/GnbhK/src/utils.jl:59</span></span>
<span class="line"><span>┌ Warning: TODO forward zero-set of arraycopy used memset rather than runtime type </span></span>
<span class="line"><span>│ Caused by:</span></span>
<span class="line"><span>│ Stacktrace:</span></span>
<span class="line"><span>│  [1] copy</span></span>
<span class="line"><span>│    @ ./array.jl:411</span></span>
<span class="line"><span>│  [2] map</span></span>
<span class="line"><span>│    @ ./tuple.jl:293</span></span>
<span class="line"><span>│  [3] map</span></span>
<span class="line"><span>│    @ ./tuple.jl:294</span></span>
<span class="line"><span>│  [4] map</span></span>
<span class="line"><span>│    @ ./namedtuple.jl:265</span></span>
<span class="line"><span>│  [5] copy</span></span>
<span class="line"><span>│    @ ~/.julia/packages/StructArrays/KLKnN/src/structarray.jl:467</span></span>
<span class="line"><span>└ @ Enzyme.Compiler ~/.julia/packages/GPUCompiler/GnbhK/src/utils.jl:59</span></span>
<span class="line"><span>┌ Warning: TODO forward zero-set of arraycopy used memset rather than runtime type </span></span>
<span class="line"><span>│ Caused by:</span></span>
<span class="line"><span>│ Stacktrace:</span></span>
<span class="line"><span>│  [1] copy</span></span>
<span class="line"><span>│    @ ./array.jl:411</span></span>
<span class="line"><span>│  [2] map</span></span>
<span class="line"><span>│    @ ./tuple.jl:293</span></span>
<span class="line"><span>│  [3] map</span></span>
<span class="line"><span>│    @ ./tuple.jl:294</span></span>
<span class="line"><span>│  [4] map</span></span>
<span class="line"><span>│    @ ./namedtuple.jl:265</span></span>
<span class="line"><span>│  [5] copy</span></span>
<span class="line"><span>│    @ ~/.julia/packages/StructArrays/KLKnN/src/structarray.jl:467</span></span>
<span class="line"><span>└ @ Enzyme.Compiler ~/.julia/packages/GPUCompiler/GnbhK/src/utils.jl:59</span></span>
<span class="line"><span>┌ Warning: TODO forward zero-set of arraycopy used memset rather than runtime type </span></span>
<span class="line"><span>│ Caused by:</span></span>
<span class="line"><span>│ Stacktrace:</span></span>
<span class="line"><span>│  [1] copy</span></span>
<span class="line"><span>│    @ ./array.jl:411</span></span>
<span class="line"><span>│  [2] map</span></span>
<span class="line"><span>│    @ ./tuple.jl:293</span></span>
<span class="line"><span>│  [3] map</span></span>
<span class="line"><span>│    @ ./tuple.jl:294</span></span>
<span class="line"><span>│  [4] map</span></span>
<span class="line"><span>│    @ ./namedtuple.jl:265</span></span>
<span class="line"><span>│  [5] copy</span></span>
<span class="line"><span>│    @ ~/.julia/packages/StructArrays/KLKnN/src/structarray.jl:467</span></span>
<span class="line"><span>└ @ Enzyme.Compiler ~/.julia/packages/GPUCompiler/GnbhK/src/utils.jl:59</span></span>
<span class="line"><span>┌ Warning: TODO forward zero-set of arraycopy used memset rather than runtime type </span></span>
<span class="line"><span>│ Caused by:</span></span>
<span class="line"><span>│ Stacktrace:</span></span>
<span class="line"><span>│  [1] copy</span></span>
<span class="line"><span>│    @ ./array.jl:411</span></span>
<span class="line"><span>│  [2] map</span></span>
<span class="line"><span>│    @ ./tuple.jl:294</span></span>
<span class="line"><span>│  [3] map</span></span>
<span class="line"><span>│    @ ./namedtuple.jl:265</span></span>
<span class="line"><span>│  [4] copy</span></span>
<span class="line"><span>│    @ ~/.julia/packages/StructArrays/KLKnN/src/structarray.jl:467</span></span>
<span class="line"><span>└ @ Enzyme.Compiler ~/.julia/packages/GPUCompiler/GnbhK/src/utils.jl:59</span></span>
<span class="line"><span>[ Info: Found initial step size 0.003125</span></span></code></pre></div><div class="warning custom-block"><p class="custom-block-title">Warning</p><p>This should be run for longer!</p></div><p>Now that we have our posterior, we can assess which parts of the image are strongly inferred by the data. This is rather unique to <code>Comrade</code> where more traditional imaging algorithms like CLEAN and RML are inherently unable to assess uncertainty in their reconstructions.</p><p>To explore our posterior let&#39;s first create images from a bunch of draws from the posterior</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">msamples </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> skymodel</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">.(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Ref</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(post), chain[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">501</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">end</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]);</span></span></code></pre></div><p>The mean image is then given by</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> StatsBase</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">imgs </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> intensitymap</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">.(msamples, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Ref</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(g))</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">mimg </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> mean</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(imgs)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">simg </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> std</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(imgs)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">fig </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> CM</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Figure</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(;resolution</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">700</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">700</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">));</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">axs </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [CM</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Axis</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(fig[i, j], xreversed</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">true</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, aspect</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">for</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> i </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">in</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, j </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">in</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">CM</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">image!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(axs[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">], mimg, colormap</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:afmhot</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">); axs[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">title</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;Mean&quot;</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">CM</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">image!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(axs[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">], simg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">./</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">max</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">.(mimg, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1e-8</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)), colorrange</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">), colormap</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:afmhot</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">);axs[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">title </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> &quot;Std&quot;</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">CM</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">image!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(axs[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">], imgs[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">],   colormap</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:afmhot</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">);</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">CM</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">image!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(axs[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">], imgs[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">end</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">], colormap</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:afmhot</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">);</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">CM</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">hidedecorations!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">.(axs)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">fig </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">|&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> DisplayAs</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">PNG </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">|&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> DisplayAs</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Text</span></span></code></pre></div><p><img src="`+t+`" alt=""></p><p>Now let&#39;s see whether our residuals look better.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">p </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Plots</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(layout</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">));</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">for</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> s </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">in</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> sample</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(chain[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">501</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">end</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">], </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">10</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">    residual!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(post, s)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">p </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">|&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> DisplayAs</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">PNG </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">|&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> DisplayAs</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Text</span></span></code></pre></div><p><img src="`+h+'" alt=""></p><p>And viola, you have a quick and preliminary image of M87 fitting only closure products. For a publication-level version we would recommend</p><ol><li><p>Running the chain longer and multiple times to properly assess things like ESS and R̂ (see <a href="/Comrade.jl/v0.11.4/tutorials/beginner/GeometricModeling#Geometric-Modeling-of-EHT-Data">Geometric Modeling of EHT Data</a>)</p></li><li><p>Fitting gains. Typically gain amplitudes are good to 10-20% for the EHT not the infinite uncertainty closures implicitly assume</p></li></ol><hr><p><em>This page was generated using <a href="https://github.com/fredrikekre/Literate.jl" target="_blank" rel="noreferrer">Literate.jl</a>.</em></p>',81)]))}const m=a(r,[["render",k]]);export{u as __pageData,m as default};
