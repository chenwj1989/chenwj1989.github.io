    
// KaTeX v0.10.2

document.addEventListener("DOMContentLoaded", function() {
    renderMathInElement(document.body, {'delimiters' : [
        {left: "$$", right: "$$", display: true},
        {left: "\\[", right: "\\]", display: true},
        {left: "$", right: "$", display: false},
        {left: "\\(", right: "\\)", display: false}
    ]});

