"""Exportação local do resumo, sem dependências externas ou chamadas de rede."""
from pathlib import Path
import base64
import html
import re
import unicodedata

ROOT = Path(__file__).resolve().parents[1]


def inline(s):
    s=html.escape(s)
    s=re.sub(r'`([^`]+)`',r'<code>\1</code>',s)
    s=re.sub(r'\*\*(.+?)\*\*',r'<strong>\1</strong>',s)
    s=re.sub(r'\[([^]]+)\]\(([^)]+)\)',r'<a href="\2">\1</a>',s)
    return s


def main():
    lines=(ROOT/'EXECUTIVE_SUMMARY.md').read_text().splitlines()
    parts=[];toc=[];i=0;in_complement=False
    while i<len(lines):
        line=lines[i].strip()
        if not line:
            i+=1;continue
        if line.startswith('#'):
            level=len(line)-len(line.lstrip('#'))
            title=line[level:].strip()
            anchor=re.sub(r'[^a-z0-9]+','-',unicodedata.normalize('NFKD',title).encode('ascii','ignore').decode().lower()).strip('-')
            if level==2:
                if in_complement:
                    parts.append('</section></details>')
                    in_complement=False
                if title.startswith('Complemento'):
                    parts.append(f'<details class="complement" id="{anchor}"><summary>{inline(title)}</summary><section>')
                    in_complement=True
                else:
                    toc.append(f'<li><a href="#{anchor}">{inline(title)}</a></li>')
                    parts.append(f'<h{level} id="{anchor}">{inline(title)}</h{level}>')
            else:
                parts.append(f'<h{level} id="{anchor}">{inline(title)}</h{level}>')
        elif line.startswith('!['):
            match=re.fullmatch(r'!\[(.*?)\]\((.*?)\)',line)
            if not match:raise ValueError(line)
            caption,path=match.groups();raw=(ROOT/path).read_bytes()
            data=base64.b64encode(raw).decode()
            parts.append(f'<figure><img alt="{html.escape(caption)}" src="data:image/png;base64,{data}"><figcaption>{html.escape(caption)}</figcaption></figure>')
        elif line.startswith('- '):
            items=[]
            while i<len(lines) and lines[i].strip().startswith('- '):
                items.append('<li>'+inline(lines[i].strip()[2:])+'</li>')
                i+=1
            parts.append('<ul>'+''.join(items)+'</ul>')
            continue
        elif line.startswith('|'):
            rows=[]
            while i<len(lines) and lines[i].strip().startswith('|'):
                cells=[x.strip() for x in lines[i].strip().strip('|').split('|')]
                if not all(re.fullmatch(r'[-: ]+',x) for x in cells):rows.append(cells)
                i+=1
            parts.append('<div class="table"><table><thead><tr>'+''.join(f'<th>{inline(x)}</th>' for x in rows[0])+'</tr></thead><tbody>'+''.join('<tr>'+''.join(f'<td>{inline(x)}</td>' for x in row)+'</tr>' for row in rows[1:])+'</tbody></table></div>')
            continue
        else:
            paragraph=[line]
            while i+1<len(lines) and lines[i+1].strip() and not lines[i+1].startswith(('#','![','|','- ')):
                i+=1;paragraph.append(lines[i].strip())
            parts.append('<p>'+inline(' '.join(paragraph))+'</p>')
        i+=1
    if in_complement:
        parts.append('</section></details>')
    style='''*{box-sizing:border-box}html{scroll-behavior:smooth}body{font:18px/1.75 system-ui,sans-serif;color:#20313b;background:#eef3f3;margin:0}main{max-width:1000px;margin:auto;padding:32px 48px 70px;background:white}h1{font-size:2.15rem;line-height:1.2;color:#126b69;margin-top:24px}h2{margin-top:3rem;padding-top:1rem;border-top:1px solid #dce7e6;color:#126b69;font-size:1.4rem;scroll-margin-top:20px}p,li{max-width:850px}strong{color:#173f48}ul{padding-left:26px}li+li{margin-top:12px}figure{margin:30px 0;background:#fff;border:1px solid #e0e8e8;border-radius:8px;padding:12px}img{display:block;width:100%;height:auto}figcaption{font-size:.87rem;color:#53616a;margin-top:12px;line-height:1.5}table{border-collapse:collapse;width:100%;font-size:.88rem;line-height:1.55}th,td{padding:12px;text-align:left;vertical-align:top;border-bottom:1px solid #d9e2e7}th{background:#e4f1f0}a{color:#126b69;text-underline-offset:3px}a:focus-visible,button:focus-visible,summary:focus-visible{outline:3px solid #bf701b;outline-offset:4px}.table{overflow-x:auto;margin:24px 0}code{background:#edf1f4;padding:2px 5px;font-size:.85em;overflow-wrap:anywhere}nav{padding:18px 24px;background:#edf6f5;border-radius:8px;margin:26px 0}nav ol{margin:10px 0;padding-left:0;list-style:none}nav li+li{margin-top:4px}nav a{font-size:.95rem}summary{cursor:pointer;font-weight:650;color:#126b69;padding:18px 0}details{border-top:1px solid #dce7e6;margin-top:24px}details section{padding-bottom:12px}.toolbar{display:flex;align-items:center;justify-content:space-between;gap:20px;color:#5f727a;font-size:.82rem}.toolbar span{text-transform:uppercase;letter-spacing:.06em}button{background:white;border:1px solid #c5d6d5;border-radius:6px;padding:9px 12px;color:#126b69;font:inherit;cursor:pointer}footer{margin-top:40px;padding-top:20px;border-top:1px solid #dce7e6;font-size:.85rem;color:#53616a}@media(max-width:650px){body{font-size:17px}main{padding:20px 18px 40px}h1{font-size:1.8rem}figure{padding:4px;margin:24px -10px}figcaption{padding:6px}th,td{padding:9px;min-width:120px}.toolbar{align-items:flex-start}nav{padding:15px}code{font-size:.8em}}@media print{@page{size:A4;margin:16mm}body{background:white;font-size:10pt}main{max-width:none;padding:0}nav,.toolbar{display:none}h1{font-size:23pt}h2{font-size:15pt;margin-top:20pt}figure{break-inside:avoid;padding:0;border:0}figure img{max-height:225mm;object-fit:contain}tr{break-inside:avoid}h2,summary{break-after:avoid}a{color:inherit}details section{display:block!important}.table{overflow:visible}footer{font-size:9pt}}'''
    navigation='<nav aria-label="Etapas da leitura"><strong>Percurso de leitura</strong><ol>'+''.join(toc)+'</ol><small>Glossário, leitura dos outros gráficos e fontes estão nos complementos ao final.</small></nav>'
    # Abertura curta antes da navegação; passos e complementos seguem na mesma ordem do Markdown.
    first_section=next(j for j,part in enumerate(parts) if part.startswith('<h2'))
    parts.insert(first_section,navigation)
    toolbar='<div class="toolbar"><span>Perfis de WPS · leitura guiada</span><button type="button" onclick="window.print()">Imprimir / salvar PDF</button></div>'
    script='''<script>let closedBeforePrint=[];window.addEventListener('beforeprint',()=>{closedBeforePrint=Array.from(document.querySelectorAll('details:not([open])'));closedBeforePrint.forEach(el=>el.open=true)});window.addEventListener('afterprint',()=>{closedBeforePrint.forEach(el=>el.open=false)})</script>'''
    document='<!doctype html><html lang="pt-BR"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Leitura guiada · Perfis de WPS</title><style>'+style+'</style></head><body><main>'+toolbar+''.join(parts)+'<footer>Documento local com figuras incorporadas. Os links de apoio usam os arquivos da pasta da análise.</footer></main>'+script+'</body></html>'
    (ROOT/'RESUMO_EXECUTIVO.html').write_text(document)
    print('Guia HTML com navegação, complementos expansíveis e imagens incorporadas.')


if __name__=='__main__':main()
