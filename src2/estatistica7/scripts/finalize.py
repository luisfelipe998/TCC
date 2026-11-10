"""Confere o pacote entregue sem recalcular experimentos ou alterar as fontes."""
from core import *
from datetime import datetime,timezone
from PIL import Image
import re,subprocess,xml.etree.ElementTree as ET


def main():
    started=json.loads((ROOT/'results/manifest_started.json').read_text())
    assert sha(ROOT/'configs/wco_v1.json')==started['config_sha256']
    for name,value in started['sources'].items():assert sha(SRC/name)==value,name
    previous=json.loads((SRC/'estatistica6/results/manifest_started.json').read_text())
    for name,value in previous['sources'].items():assert sha(SRC/name)==value,name
    assert json.loads((ROOT/'validation/final.json').read_text())['status']=='PASS'
    assert json.loads((ROOT/'results/phase_null_metadata.json').read_text())['replicates_completed']==199
    pngs=sorted((ROOT/'figures').glob('*.png'));assert len(pngs)==9
    for path in pngs:
        with Image.open(path) as im:
            assert min(im.size)>=500;im.verify()
        assert path.with_suffix('.pdf').stat().st_size>1000
    links=[]
    for path in ROOT.glob('*.md'):
        for target in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)',path.read_text()):
            if '://' in target or target.startswith('#'):continue
            local=(path.parent/target.split('#')[0]).resolve()
            if local!=ROOT/'results/artifact_manifest.json':assert local.exists(),(path.name,target)
            links.append(dict(document=path.name,target=target))
    pdf=ROOT/'RESUMO_EXECUTIVO.pdf'
    xml=subprocess.check_output(['pdftotext','-bbox-layout',str(pdf),'-'])
    tree=ET.fromstring(xml);pages=[node for node in tree.iter() if node.tag.endswith('}page')]
    assert len(pages)==4
    word_count=0
    for page in pages:
        width,height=float(page.attrib['width']),float(page.attrib['height'])
        for word in page.iter():
            if not word.tag.endswith('}word'):continue
            word_count+=1
            assert 0<=float(word.attrib['xMin'])<=float(word.attrib['xMax'])<=width+.1
            assert 0<=float(word.attrib['yMin'])<=float(word.attrib['yMax'])<=height+.1
    content=subprocess.check_output(['pdftotext','-layout',str(pdf),'-']).decode()
    for phrase in ['0,986','0,979','68,3%','22 representações','0,69','47,2%','Plug e slug']:
        assert phrase in content,phrase
    # Conferência da confusão operacional citada nos documentos.
    inventory=csv_read(ROOT/'results/inventory.csv')
    assert (inventory.groupby(['m_code','q_code']).pattern.nunique()==1).all()
    q=csv_read(ROOT/'results/within_q.csv')
    estimable=q[q.estimable][['scope','q_code','n_a','n_b']].drop_duplicates()
    assert len(estimable)==1 and estimable.iloc[0].scope=='churn_annular' and estimable.iloc[0].q_code==40
    report=dict(status='PASS',pdf_pages=len(pages),pdf_text_words=word_count,pdf_text_inside_page=True,
                png_pdf_figure_pairs=len(pngs),local_markdown_links_checked=len(links),
                direct_sources_unchanged=len(started['sources']),prior_source_hashes_verified=len(previous['sources']),
                scope='Integridade dos arquivos e limites do texto no PDF; a inspeção visual e a auditoria numérica são complementares.')
    json_write(ROOT/'validation/artifacts.json',report)
    manifest_path=ROOT/'results/artifact_manifest.json'
    paths=sorted(p for p in ROOT.rglob('*') if p.is_file() and p!=manifest_path and '__pycache__' not in p.parts and p.name!='.DS_Store')
    manifest=dict(completed_utc=datetime.now(timezone.utc).isoformat(),source_root=str(SRC),sources=started['sources'],
                  artifacts={str(p.relative_to(ROOT)):sha(p) for p in paths},
                  note='O próprio manifesto é excluído da lista de hashes; arquivos de cache Python não são artefatos da análise.')
    json_write(manifest_path,manifest)
    assert all(sha(ROOT/name)==digest for name,digest in manifest['artifacts'].items())
    print(json.dumps(dict(**report,artifact_hashes_verified=len(paths)),ensure_ascii=False,indent=2),flush=True)


if __name__=='__main__':main()
