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


# Tutorials

## Beginner Tutorials

<Gallery :images="beginner" />
```