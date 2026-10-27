"""Chemical formula parsing, standing in for the StoichTools parse_formula
used by setupphases.m.

``parse_formula("Ca2(Mg0.75Fe0.25)5Si8O22(OH)2")`` returns element counts in
order of first appearance: ``{"Ca": 2, "Mg": 3.75, "Fe": 1.25, "Si": 8, "O": 24, "H": 2}``.
"""

import re

_TOKEN = re.compile(r"([A-Z][a-z]?)|(\()|(\))|(\d+(?:\.\d*)?|\.\d+)")


def parse_formula(formula):
    """Element counts for a formula with nested parentheses and decimal subscripts."""
    s = re.sub(r"\s", "", formula)
    stack = [{}]
    pos = 0
    last = None  # (kind, value) of the item a following number multiplies

    def add(d, el, n):
        d[el] = d.get(el, 0.0) + n

    while pos < len(s):
        m = _TOKEN.match(s, pos)
        if not m:
            raise ValueError(f"cannot parse formula {formula!r} at {s[pos:]!r}")
        pos = m.end()
        el, op, cl, num = m.groups()
        if el:
            add(stack[-1], el, 1.0)
            last = ("el", el)
        elif op:
            stack.append({})
            last = None
        elif cl:
            group = stack.pop()
            for k, v in group.items():
                add(stack[-1], k, v)
            last = ("group", group)
        else:
            n = float(num)
            if last is None:
                raise ValueError(f"number without element in {formula!r}")
            if last[0] == "el":
                add(stack[-1], last[1], n - 1)
            else:
                for k, v in last[1].items():
                    add(stack[-1], k, v * (n - 1))
            last = None
    if len(stack) != 1:
        raise ValueError(f"unbalanced parentheses in {formula!r}")
    return stack[0]


def hydrate_to_formula(formula):
    """setstoichiometry's handling of ``X:nH2O``: returns ``X(H2O)n`` without spaces."""
    f, _, rest = formula.partition(":")
    rest = rest.split(":")[0] if rest else ""
    if rest:
        i = rest.find("H")
        s, h2o = (rest[:i], rest[i:]) if i >= 0 else (rest, "")
        if s:
            f = f"{f}({h2o}){s}"
    return re.sub(r"\s", "", f)
