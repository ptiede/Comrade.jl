// .vitepress/theme/index.ts
import { h } from 'vue'
import type { Theme } from 'vitepress'
import DefaultTheme from 'vitepress/theme'
import AsideRelated from '../../components/AsideRelated.vue'

import { enhanceAppWithTabs } from 'vitepress-plugin-tabs/client'
import './style.css'

export default {
  extends: DefaultTheme,
  Layout() {
    return h(DefaultTheme.Layout, null, {
      // 'home-hero-info-after': () => h(HomeTrustees),
      'aside-ads-before': () => h(AsideRelated),
    })
  },
  enhanceApp({ app, router, siteData }) {
    enhanceAppWithTabs(app)
  }
} satisfies Theme