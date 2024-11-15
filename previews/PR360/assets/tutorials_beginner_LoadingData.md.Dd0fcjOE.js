import{_ as s,c as i,o as a,a6 as t}from"./chunks/framework.D_Wlhr2D.js";const n="/Comrade.jl/previews/PR360/assets/yavysfl.Dn8CuuDF.png",e="/Comrade.jl/previews/PR360/assets/egrhifn.DTHJMRBG.png",l="/Comrade.jl/previews/PR360/assets/tolgznh.IQgfKv3m.png",F=JSON.parse('{"title":"Loading Data into Comrade","description":"","frontmatter":{},"headers":[],"relativePath":"tutorials/beginner/LoadingData.md","filePath":"tutorials/beginner/LoadingData.md","lastUpdated":null}'),p={name:"tutorials/beginner/LoadingData.md"},h=t(`<h1 id="Loading-Data-into-Comrade" tabindex="-1">Loading Data into Comrade <a class="header-anchor" href="#Loading-Data-into-Comrade" aria-label="Permalink to &quot;Loading Data into Comrade {#Loading-Data-into-Comrade}&quot;">​</a></h1><p>The VLBI field does not have a standardized data format, and the EHT uses a particular uvfits format similar to the optical interferometry oifits format. As a result, we reuse the excellent <code>eht-imaging</code> package to load data into <code>Comrade</code>.</p><p>Once the data is loaded, we then convert the data into the tabular format <code>Comrade</code> expects. Note that this may change to a Julia package as the Julia radio astronomy group grows.</p><p>To get started, we will load <code>Comrade</code> and <code>Plots</code> to enable visualizations of the data</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Comrade</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Plots</span></span></code></pre></div><p>We also load Pyehtim since it loads eht-imaging into Julia using PythonCall and exports the variable ehtim</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Pyehtim</span></span></code></pre></div><p>To load the data we will use <code>eht-imaging</code>. We will use the 2017 public M87 data which can be downloaded from <a href="https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/EHTC_FirstM87Results_Apr2019" target="_blank" rel="noreferrer">cyverse</a></p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">obseht </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> ehtim</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">obsdata</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">load_uvfits</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">joinpath</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(__DIR, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;..&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;..&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;Data&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;SR1_M87_2017_096_lo_hops_netcal_StokesI.uvfits&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">))</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>Python: &lt;ehtim.obsdata.Obsdata object at 0x7f4960b40820&gt;</span></span></code></pre></div><p>Now we will average the data over telescope scans. Note that the EHT data has been pre-calibrated so this averaging doesn&#39;t induce large coherence losses.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">obs </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Pyehtim</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">scan_average</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(obseht)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>Python: &lt;ehtim.obsdata.Obsdata object at 0x7f497b3156f0&gt;</span></span></code></pre></div><div class="warning custom-block"><p class="custom-block-title">Warning</p><p>We use a custom scan-averaging function to ensure that the scan-times are homogenized.</p></div><p>We can now extract data products that <code>Comrade</code> can use</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">vis    </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> extract_table</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(obs, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Visibilities</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">()) </span><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">## complex visibilites</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">amp    </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> extract_table</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(obs, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">VisibilityAmplitudes</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">()) </span><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">## visibility amplitudes</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">cphase </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> extract_table</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(obs, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">ClosurePhases</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(; snrcut</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">3.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)) </span><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">## extract minimal set of closure phases</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">lcamp  </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> extract_table</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(obs, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">LogClosureAmplitudes</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(; snrcut</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">3.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)) </span><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">## extract minimal set of log-closure amplitudes</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>EHTObservationTable{Comrade.EHTLogClosureAmplitudeDatum{:I}}</span></span>
<span class="line"><span>  source:      M87</span></span>
<span class="line"><span>  mjd:         57849</span></span>
<span class="line"><span>  bandwidth:   1.856e9</span></span>
<span class="line"><span>  sites:       [:AA, :AP, :AZ, :JC, :LM, :PV, :SM]</span></span>
<span class="line"><span>  nsamples:    129</span></span></code></pre></div><p>For polarization we first load the data in the cirular polarization basis Additionally, we load the array table at the same time to load the telescope mounts.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">obseht </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Pyehtim</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">load_uvfits_and_array</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">                joinpath</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(__DIR, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;..&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;..&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;Data&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;polarized_gaussian_all_corruptions.uvfits&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">),</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">                joinpath</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(__DIR, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;..&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;..&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;Data&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;array.txt&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">),</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">                polrep</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;circ&quot;</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">                        )</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">obs </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Pyehtim</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">scan_average</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(obseht)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">coh </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> extract_table</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(obs, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Coherencies</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">())</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>EHTObservationTable{Comrade.EHTCoherencyDatum{Float64}}</span></span>
<span class="line"><span>  source:      17.761122472222223:-28.992189444444445</span></span>
<span class="line"><span>  mjd:         51544</span></span>
<span class="line"><span>  bandwidth:   1.0e9</span></span>
<span class="line"><span>  sites:       [:AA, :AP, :AZ, :JC, :LM, :PV, :SM]</span></span>
<span class="line"><span>  nsamples:    315</span></span></code></pre></div><div class="warning custom-block"><p class="custom-block-title">Warning</p><p>Always use our <code>extract_cphase</code> and <code>extract_lcamp</code> functions to find the closures eht-imaging will sometimes incorrectly calculate a non-redundant set of closures.</p></div><p>We can also recover the array used in the observation using</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> DisplayAs</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ac </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> arrayconfig</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(vis)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(ac) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">|&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> DisplayAs</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">PNG </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">|&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> DisplayAs</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Text </span><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># Plot the baseline coverage</span></span></code></pre></div><p><img src="`+n+`" alt=""></p><p>To plot the data we just call</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">l </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> @layout</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [a b; c d]</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">pv </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> plot</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(vis);</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">pa </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> plot</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(amp);</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">pcp </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> plot</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(cphase);</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">plc </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> plot</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(lcamp);</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(pv, pa, pcp, plc; layout</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">l) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">|&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> DisplayAs</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">PNG </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">|&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> DisplayAs</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Text</span></span></code></pre></div><p><img src="`+e+'" alt=""></p><p>And also the coherency matrices</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(coh) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">|&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> DisplayAs</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">PNG </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">|&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> DisplayAs</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Text</span></span></code></pre></div><p><img src="'+l+'" alt=""></p><hr><p><em>This page was generated using <a href="https://github.com/fredrikekre/Literate.jl" target="_blank" rel="noreferrer">Literate.jl</a>.</em></p>',32),k=[h];function d(r,o,E,g,c,y){return a(),i("div",null,k)}const b=s(p,[["render",d]]);export{F as __pageData,b as default};