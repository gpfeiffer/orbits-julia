-- Pandoc filter: <div class="alert alert-KIND"> ... </div>  ->  nbalert box.
--
-- Pandoc's native_divs extension already parses the notebooks' alert divs into
-- Div elements (classes "alert", "alert-danger", possibly "alert-block"); the
-- LaTeX writer would just drop the div.  Wrap its content in the nbalert
-- environment that templates/book/index.tex.j2 defines instead.
local kinds = { danger = true, warning = true, info = true, success = true }

function Div(el)
  for _, class in ipairs(el.classes) do
    local kind = class:match('^alert%-(%a+)$')
    if kind and kinds[kind] then
      return {
        pandoc.RawBlock('latex', '\\begin{nbalert}{' .. kind .. '}'),
        el,
        pandoc.RawBlock('latex', '\\end{nbalert}'),
      }
    end
  end
end
