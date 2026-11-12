"""Synthetic WebBook-like pages for the tests: the WebBook's markup, made-up
values and references (no NIST data is redistributed)."""


def page(title, formula_html, cas, rows, refs):
    """rows: (class, quantity, value, units, method, ref_id, short, comment)"""
    trs = "\n".join(
        f'<tr class="{c}"><td>{q}</td><td>{v}</td><td>{u}</td><td>{m}</td>'
        f'<td><a href="#{rid}">{short}</a></td><td>{comment}</td></tr>'
        for c, q, v, u, m, rid, short, comment in rows
    )
    spans = "\n".join(
        f'<span id="{rid}">{text} [ all data ]</span><br>' for rid, text in refs
    )
    return f"""<html><head><title>{title}</title></head><body>
<ul><li><strong><a href="#">Formula</a>:</strong> {formula_html}</li>
<li><strong>CAS Registry Number:</strong> {cas}</li></ul>
<h2 id="Thermo-Gas">Gas phase thermochemistry data</h2>
<table>{trs}</table>
<hr>
<h2 id="Refs">References</h2>
<p>{spans}</p>
</body></html>"""
