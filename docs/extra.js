/* pkgdown/extra.js
 * Replaces pkgdown's built-in floating "On this page" TOC (powered by
 * bootstrap-toc, which only supports ~2 nested levels and shows/hides
 * nested items abruptly while scrolling) with a custom TOC that:
 *   - shows ALL heading levels (h2-h6) found on the page
 *   - lets the user click to expand/collapse any branch
 *   - only highlights the current section while scrolling, instead of
 *     hiding and re-showing branches
 *
 * Just drop this file at pkgdown/extra.js in your package root and
 * rebuild the site (pkgdown::build_site()). No _pkgdown.yml changes
 * are required, though setting `toc: depth: 5` there too doesn't hurt.
 */
(function () {
  "use strict";

  function ready(fn) {
    if (document.readyState !== "loading") fn();
    else document.addEventListener("DOMContentLoaded", fn);
  }

  ready(function () {
    // pkgdown marks its generated floating TOC with data-toggle="toc"
    // (this attribute name has stayed stable across BS3/4/5 templates).
    var tocNav = document.querySelector('nav[data-toggle="toc"]');
    if (!tocNav) return;

    // Main article content: try the most specific pkgdown container
    // first, then fall back to broader selectors.
    var content =
      document.querySelector("#main .contents") ||
      document.querySelector("#main") ||
      document.querySelector("main .contents") ||
      document.querySelector("main");
    if (!content) return;

    var headings = content.querySelectorAll(
      "h2[id]:not([data-toc-skip]), h3[id]:not([data-toc-skip]), " +
        "h4[id]:not([data-toc-skip]), h5[id]:not([data-toc-skip]), " +
        "h6[id]:not([data-toc-skip])"
    );
    if (!headings.length) return;

    // --- Stop pkgdown/bootstrap-toc from touching this nav any more ---
    tocNav.removeAttribute("data-toggle");
    if (window.bootstrap && window.bootstrap.ScrollSpy) {
      var instance = window.bootstrap.ScrollSpy.getInstance(document.body);
      if (instance) instance.dispose();
    }

    // --- Build a nested tree from the flat heading list ---
    var root = { level: 1, children: [] };
    var stack = [root];
    headings.forEach(function (h) {
      var level = parseInt(h.tagName.substring(1), 10);
      var node = {
        level: level,
        id: h.id,
        text: h.textContent.trim(),
        children: [],
      };
      while (stack.length > 1 && stack[stack.length - 1].level >= level) {
        stack.pop();
      }
      stack[stack.length - 1].children.push(node);
      stack.push(node);
    });

    function buildList(nodes) {
      var ul = document.createElement("ul");
      nodes.forEach(function (node) {
        var li = document.createElement("li");
        var hasChildren = node.children.length > 0;
        li.className = hasChildren ? "vtoc-parent" : "vtoc-leaf";

        var toggle = document.createElement("span");
        toggle.className = "vtoc-toggle";
        toggle.setAttribute("aria-hidden", "true");
        li.appendChild(toggle);

        var a = document.createElement("a");
        a.href = "#" + node.id;
        a.textContent = node.text;
        li.appendChild(a);

        if (hasChildren) {
          var childUl = buildList(node.children);
          childUl.className = "vtoc-collapsed";
          li.appendChild(childUl);

          function toggleBranch(e) {
            e.preventDefault();
            var open = li.classList.toggle("vtoc-open");
            childUl.classList.toggle("vtoc-collapsed", !open);
          }
          toggle.addEventListener("click", toggleBranch);
          toggle.tabIndex = 0;
          toggle.setAttribute("role", "button");
          toggle.addEventListener("keydown", function (e) {
            if (e.key === "Enter" || e.key === " ") toggleBranch(e);
          });
        }
        ul.appendChild(li);
      });
      return ul;
    }

    var wrapper = document.createElement("div");
    wrapper.className = "vtoc";

    var controls = document.createElement("div");
    controls.className = "vtoc-controls";
    var expandAll = document.createElement("button");
    expandAll.type = "button";
    expandAll.textContent = "Expand all";
    var collapseAll = document.createElement("button");
    collapseAll.type = "button";
    collapseAll.textContent = "Collapse all";
    controls.appendChild(expandAll);
    controls.appendChild(collapseAll);

    var list = buildList(root.children);

    wrapper.appendChild(controls);
    wrapper.appendChild(list);

    // Replace whatever pkgdown put inside the nav (keep the "On this
    // page" heading if pkgdown already rendered one just before it).
    tocNav.innerHTML = "";
    var heading = document.createElement("div");
    heading.className = "vtoc-heading";
    heading.textContent = "On this page";
    tocNav.appendChild(heading);
    tocNav.appendChild(wrapper);

    expandAll.addEventListener("click", function () {
      wrapper.querySelectorAll("li.vtoc-parent").forEach(function (li) {
        li.classList.add("vtoc-open");
      });
      wrapper.querySelectorAll("ul.vtoc-collapsed").forEach(function (ul) {
        ul.classList.remove("vtoc-collapsed");
      });
    });
    collapseAll.addEventListener("click", function () {
      wrapper.querySelectorAll("li.vtoc-parent").forEach(function (li) {
        li.classList.remove("vtoc-open");
      });
      wrapper.querySelectorAll("ul").forEach(function (ul) {
        if (ul !== list) ul.classList.add("vtoc-collapsed");
      });
    });

    // --- Highlight current section on scroll, without hiding others ---
    var linkFor = {};
    wrapper.querySelectorAll("a[href^='#']").forEach(function (a) {
      linkFor[a.getAttribute("href").substring(1)] = a;
    });

    function setActive(id) {
      wrapper.querySelectorAll("a.vtoc-active").forEach(function (a) {
        a.classList.remove("vtoc-active");
      });
      var a = linkFor[id];
      if (!a) return;
      a.classList.add("vtoc-active");
      // make sure every ancestor branch is expanded so the active
      // item is actually visible
      var li = a.closest("li");
      while (li) {
        var parentUl = li.parentElement;
        if (parentUl && parentUl.classList.contains("vtoc-collapsed")) {
          parentUl.classList.remove("vtoc-collapsed");
          var parentLi = parentUl.closest("li");
          if (parentLi) parentLi.classList.add("vtoc-open");
        }
        li = li.parentElement ? li.parentElement.closest("li") : null;
      }
    }

    if ("IntersectionObserver" in window) {
      var observer = new IntersectionObserver(
        function (entries) {
          var visible = entries.filter(function (e) {
            return e.isIntersecting;
          });
          if (visible.length) {
            visible.sort(function (a, b) {
              return a.boundingClientRect.top - b.boundingClientRect.top;
            });
            setActive(visible[0].target.id);
          }
        },
        { rootMargin: "0px 0px -70% 0px", threshold: [0, 1] }
      );
      headings.forEach(function (h) {
        observer.observe(h);
      });
    }
  });
})();
