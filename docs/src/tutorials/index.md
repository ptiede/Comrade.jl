```@raw html
<script setup lang="ts">
import Gallery from "../components/Gallery.vue";

const beginner = [
  {
    href: "beginner/LoadingData", 
    caption: "Loading Data with Pyehtim", 
    desc: "How to load data using standard eht-imaging in Julia." 
  },
  {
      href: "beginner/GeometricModeling.md",
      caption: "Geometric Modeling of M87",
      desc: "Fitting simple geometric models to the 2017 M87 data."
  }
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
</script>



# Tutorials

## Beginner Tutorials

<Gallery :images="beginner" />
```
