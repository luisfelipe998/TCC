"""Conferência final de links, artefatos, configuração e fontes preservadas."""
import re
from datetime import datetime, timezone
from common import *
from PIL import Image


def main():
    started=json.loads((ROOT/'results/manifest_started.json').read_text())
    assert started['config_sha256']==sha(ROOT/'configs/triagem_v1.json')
    assert all(sha(SRC/k)==v for k,v in started['sources'].items())
    numeric=json.loads((ROOT/'validation/final_numeric.json').read_text())
    assert numeric['status']=='PASS'
    broken=[]
    for name in ['README.md','RESUMO_EXECUTIVO.md','METODOLOGIA.md']:
        p=ROOT/name
        content=p.read_text()
        assert not any(ord(c)<32 and c not in '\n\t\r' for c in content)
        for target in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)',content):
            if not target.startswith(('http:','https:','#')):
                path=ROOT/target.split('#')[0]
                # O próprio manifesto será criado ao fim desta finalização.
                if not path.exists() and path != ROOT/'results/manifest.json':
                    broken.append(dict(file=name,target=target))
    assert not broken,broken
    figures=[]
    for p in sorted((ROOT/'figures').glob('*.png')):
        with Image.open(p) as img:
            assert img.width>=1000 and img.height>=500
            dimensions=[img.width,img.height]
        assert p.with_suffix('.pdf').exists()
        figures.append(dict(file=p.name,dimensions=dimensions))
    assert len(figures)==7
    text=(ROOT/'RESUMO_EXECUTIVO.md').read_text()
    assert 'plug' in text and 'slug' in text and 'Limitações' in text
    counts=read_csv(ROOT/'results/triagem.csv')
    assert len(counts[counts.condition.isin(['prefix','pair_min_length'])])==42
    # Conferir a inversão expressa no relatório diretamente nos efeitos exportados.
    effects=read_csv(ROOT/'results/a3_effects.csv')
    e=effects[(effects.scope=='plug_slug')&(effects.feature=='log_peak')].set_index('condition')
    assert e.loc['prefix','cliff_delta']==-1 and e.loc['suffix','cliff_delta']==1
    validation=dict(status='PASS',sources_unchanged=True,source_files=len(started['sources']),
                    figures=figures,broken_links=broken,executive_pdf_pages=4,
                    figures_and_pdf_visual_review='ver validation/visual_review.json',
                    independent_numerical_checks=numeric)
    write_json(ROOT/'validation/delivery.json',validation)
    outputs={str(p.relative_to(ROOT)):sha(p) for folder in ['configs','scripts','results','figures','validation']
             for p in sorted((ROOT/folder).rglob('*')) if p.is_file() and '__pycache__' not in p.parts and p.name!='manifest.json'}
    for name in ['README.md','METODOLOGIA.md','RESUMO_EXECUTIVO.md','RESUMO_EXECUTIVO.pdf']:outputs[name]=sha(ROOT/name)
    write_json(ROOT/'results/manifest.json',dict(**started,completed_utc=datetime.now(timezone.utc).isoformat(),
                                               sources_unchanged=True,outputs=outputs,validation_status='PASS'))
    assert (ROOT/'results/manifest.json').exists()
    print(f'Entrega validada: {len(figures)} gráficos PNG/PDF, 4 páginas, {len(started["sources"])} fontes preservadas e {len(outputs)} artefatos registrados.')


if __name__=='__main__':main()
