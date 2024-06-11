```@raw html
<script setup lang="ts">
import Gallery from "../components/Gallery.vue";

const beginner = [
  {
    href: "beginner/LoadingData", 
    src: "../assets/vis.png",
    caption: "Loading Data with Pyehtim",
    desc: "How to load data using standard eht-imaging in Julia." 
  },
  {
    href: "beginner/GeometricModeling", 
    src: "../assets/geom_model.png",
    caption: "Geometric Modeling of M87*",
    desc: "Modeling a black hole with simple geometric models" 
  }
];

</script>



# Tutorials

## Beginner Tutorials

<Gallery :images="beginner" />


```
