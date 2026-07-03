(function () {
  const storageKey = "econpdes-docs-menu-scroll";
  const maxAgeMs = 15000;

  function pageKey(url) {
    return url.pathname + url.search;
  }

  function docsMenu() {
    return document.querySelector("#documenter .docs-menu");
  }

  function activeTocItem() {
    const active = document.querySelector("#documenter .docs-menu li.is-active");
    if (!active) return null;

    for (const child of active.children) {
      if (child.classList && child.classList.contains("tocitem")) return child;
    }
    return active;
  }

  function savedScrollTop() {
    try {
      const raw = window.sessionStorage.getItem(storageKey);
      if (!raw) return null;

      const saved = JSON.parse(raw);
      if (saved.path !== pageKey(window.location)) return null;
      if (Date.now() - saved.time > maxAgeMs) return null;

      const scrollTop = Number(saved.scrollTop);
      return Number.isFinite(scrollTop) ? scrollTop : null;
    } catch (_error) {
      return null;
    }
  }

  function restoreSidebarScroll() {
    const menu = docsMenu();
    if (!menu) return;

    const saved = savedScrollTop();
    if (saved !== null) {
      menu.scrollTop = saved;
      return;
    }

    const active = activeTocItem();
    if (!active) return;

    const menuRect = menu.getBoundingClientRect();
    const activeRect = active.getBoundingClientRect();
    const activeBottom = menu.scrollTop + activeRect.bottom - menuRect.top;
    const margin = 16;

    menu.scrollTop = Math.max(0, activeBottom - menu.clientHeight + margin);
  }

  function installScrollPreserver() {
    const menu = docsMenu();
    if (!menu) return;

    menu.addEventListener("click", function (event) {
      const link = event.target.closest("a.tocitem[href]");
      if (!link) return;

      const target = new URL(link.getAttribute("href"), window.location.href);
      if (target.protocol !== window.location.protocol || target.host !== window.location.host) return;
      if (pageKey(target) === pageKey(window.location)) return;

      window.sessionStorage.setItem(
        storageKey,
        JSON.stringify({
          path: pageKey(target),
          scrollTop: menu.scrollTop,
          time: Date.now(),
        }),
      );
    });
  }

  function scheduleRestore() {
    restoreSidebarScroll();
    setTimeout(restoreSidebarScroll, 50);
    setTimeout(restoreSidebarScroll, 200);
  }

  document.addEventListener("DOMContentLoaded", function () {
    installScrollPreserver();
    scheduleRestore();
  });
  window.addEventListener("load", scheduleRestore);
})();
