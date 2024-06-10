```@raw html
<script setup lang="ts">
import Gallery from "../components/Gallery.vue";

const beginner = [
    {
        href: "beginner/LoadingData.md",
        caption: "Loading Data with Pyehtim",
        src: "../assets/vis.png",
        desc: "How to load data using standard eht-imaging in Julia."
    },
];

const intermediate = [
    {
        href: "intermediate/ClosureImaging.md",
        caption: "Imaging a Black Hole using Closures",
        desc: "Learn how to image M87* using only closure quantities."
    },
    {
        href: "intermediate/StokesIImaging.md",
        caption: "Stokes I Simultaneous Image and Instrument Modeling",
        desc: "Learn how to image complex visibilities using Comrade"
    },
    {
        href: "intermediate/PolarizedImaging",
        caption: "Full Stokes Imaging and Calibration using the RIME formalism",
        desc: "Learn how to do solve the RIME with Comrade"
    }
];

const advanced = [
    {
        href: "advanced/HybridImaging.md",
        caption: "Simultaneous modeling with residual imaging for M87*"
        desc: "Learn how to using a hybrid model decomposition of a black hole for improved inference"
    }
];

# Tutorials

## Beginner Tutorials

<Gallery :images="beginner" />

## Intermediate Tutorials

<Gallery :images="intermediate" />

## Advanced Tutorials

<Gallery :images="advanced" />

```