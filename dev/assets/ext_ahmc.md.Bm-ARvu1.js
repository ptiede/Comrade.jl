import{_ as i,o as a,c as e,a5 as t}from"./chunks/framework.kib6YrKW.js";const c=JSON.parse('{"title":"AdvancedHMC Extension","description":"","frontmatter":{},"headers":[],"relativePath":"ext/ahmc.md","filePath":"ext/ahmc.md","lastUpdated":null}'),n={name:"ext/ahmc.md"};function l(h,s,p,d,o,r){return a(),e("div",null,s[0]||(s[0]=[t(`<h1 id="AdvancedHMC-Extension" tabindex="-1">AdvancedHMC Extension <a class="header-anchor" href="#AdvancedHMC-Extension" aria-label="Permalink to &quot;AdvancedHMC Extension {#AdvancedHMC-Extension}&quot;">​</a></h1><p>The first choice when sampling from the model/image posterior, is <a href="https://github.com/TuringLang/AdvancedHMC.jl" target="_blank" rel="noreferrer"><code>AdvancedHMC</code></a>, which uses Hamiltonian Monte Carlo to sample from the posterior. Specifically, we usually use the <code>NUTS</code> algorithm.</p><p>The interface to <code>AdvancedHMC</code> is very powerful and general. To simplify the procedure for <code>Comrade</code> users, we have provided a thin interface. A user needs to specify a <code>sampler</code> and then call the <code>sample</code> function.</p><p>To sample a user can use follow the standard <code>AdvancedHMC</code> interface, e.g.,</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">chain </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> sample</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(post, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">NUTS</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.8</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">), </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">10_000</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">; n_adapts</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">5_000</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><div class="warning custom-block"><p class="custom-block-title">Warning</p><p>To use HMC the <code>VLBIPosterior</code> must be created with a specific <code>admode</code> specified. The <code>admode</code> can be a union of <code>Nothing</code> and <code>&lt;:EnzymeCore.Mode</code> types. We recommend using <code>Enzyme.set_runtime_activity(Enzyme.Reverse)</code></p></div><p>In addition our sample call has a few additional keyword arguments:</p><ul><li><code>saveto = MemoryStore()</code>: Specifies how to store the samples. The default is <code>MemoryStore</code> which stores the samples directly in RAM. For large models this is not a good idea. To save samples periodically to disk use <a href="/Comrade.jl/dev/api#Comrade.DiskStore"><code>DiskStore</code></a>, and then load the results with <code>load_samples</code>.</li></ul><p>Note that like most <code>AbstractMCMC</code> samplers the initial location can be specified with the <code>initial_params</code> argument.</p><h2 id="example" tabindex="-1">Example <a class="header-anchor" href="#example" aria-label="Permalink to &quot;Example&quot;">​</a></h2><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Comrade</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> AdvancedHMC</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Enzyme</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># Some stuff to create a posterior object</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">post </span><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># of type Comrade.Posterior</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">out </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> sample</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(post, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">NUTS</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.9</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">), </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2_000</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">; n_adapts</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1_000</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, saveto</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">DiskStore</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">())</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">chain </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> load_samples</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(out)</span></span></code></pre></div>`,11)]))}const E=i(n,[["render",l]]);export{c as __pageData,E as default};