#!/usr/bin/env python3
"""Keep the familiar v15.1 report paths on the latest verified investigation."""
from pathlib import Path
from urllib.parse import urlsplit, unquote, quote
import os
import re

BASE=Path(__file__).resolve().parents[1]
LATEST=BASE/'all_scenarios_adversarial'
NAMES=('EXECUTIVE_SUMMARY.md','EXECUTIVE_SUMMARY_PT_BR.md','FIGURE_CAPTIONS_PT_BR.md')
LINK=re.compile(r'(\[.*?\]\()(.+?)(\))')


def archive_name(name):
    return name.removesuffix('.md')+'_2026-10-23.md'


def latest_available():
    return all((LATEST/'findings'/name).is_file() for name in NAMES)


def rebase(text,source,destination,targets=None):
    targets=targets or {}
    def replace(match):
        target=match[2]
        url=urlsplit(target)
        if not url.path or url.scheme or url.netloc or target.startswith('/'):
            return match[0]
        actual=(source/unquote(url.path)).resolve()
        actual=targets.get(actual,actual)
        new=quote(Path(os.path.relpath(actual,destination)).as_posix(),safe='/._-')
        if url.query:new+='?'+url.query
        if url.fragment:new+='#'+url.fragment
        return match[1]+new+match[3]
    return LINK.sub(replace,text)


def write_legacy_report(name,lines):
    """Rebuild the earlier analysis without replacing the current entry point."""
    text=chr(10).join(lines)+chr(10)
    folder=BASE/'findings'
    if latest_available():
        targets={(folder/n).resolve():(folder/archive_name(n)).resolve() for n in NAMES}
        text=rebase(text,folder,folder,targets)
        (folder/archive_name(name)).write_text(text)
    else:
        (folder/name).write_text(text)


def latest_text(name):
    folder=BASE/'findings'
    targets={(LATEST/'findings'/n).resolve():(folder/n).resolve() for n in NAMES}
    source=LATEST/'findings'/name
    text=rebase(source.read_text(),source.parent,folder,targets)
    first,rest=text.split(chr(10),1)
    note=(f'**Atualizado em 24/10/2026 com a investigação adversarial.** [Versão anterior de 23/10](./{archive_name(name)}).'
          if name.endswith('_PT_BR.md') else
          f'**Updated on 24 October 2026 with the adversarial investigation.** [Previous version from 23 October](./{archive_name(name)}).')
    return first+chr(10)*2+note+chr(10)+rest


def check_publication():
    if not latest_available():return []
    checked=[]
    for name in NAMES:
        path=BASE/'findings'/name
        assert path.read_text()==latest_text(name),('Stale main report',path)
        assert (path.parent/archive_name(name)).exists()
        checked.append(str(path.relative_to(BASE)))
    return checked


def publish_latest():
    if not latest_available():return []
    folder=BASE/'findings'
    updated=[]
    for name in NAMES:
        archive=folder/archive_name(name)
        if not archive.exists():
            old=folder/name
            assert old.exists(),old
            assert 'all_26_selected_wall_' not in old.read_text(),('Missing original archive',old)
            old_targets={(folder/n).resolve():(folder/archive_name(n)).resolve() for n in NAMES}
            archive.write_text(rebase(old.read_text(),folder,folder,old_targets))
        (folder/name).write_text(latest_text(name))
        updated.append(str((folder/name).relative_to(BASE)))
    return updated


if __name__=='__main__':
    print(chr(10).join(publish_latest()))
