import ast
import io
import sys
import tokenize
from pathlib import Path

def remove_docstrings(source: str) -> str:
    tree = ast.parse(source)

    class DocstringRemover(ast.NodeTransformer):

        def _strip(self, node):
            if node.body and isinstance(node.body[0], ast.Expr) and isinstance(node.body[0].value, ast.Constant) and isinstance(node.body[0].value.value, str):
                node.body.pop(0)
                if not node.body:
                    node.body.append(ast.Pass())
            return node

        def visit_Module(self, node):
            self.generic_visit(node)
            if node.body and isinstance(node.body[0], ast.Expr) and isinstance(node.body[0].value, ast.Constant) and isinstance(node.body[0].value.value, str):
                node.body.pop(0)
            return node

        def visit_FunctionDef(self, node):
            self.generic_visit(node)
            return self._strip(node)

        def visit_AsyncFunctionDef(self, node):
            self.generic_visit(node)
            return self._strip(node)

        def visit_ClassDef(self, node):
            self.generic_visit(node)
            return self._strip(node)
    tree = DocstringRemover().visit(tree)
    ast.fix_missing_locations(tree)
    return ast.unparse(tree)

def remove_comments(source: str) -> str:
    out_tokens = []
    tokens = tokenize.generate_tokens(io.StringIO(source).readline)
    for tok in tokens:
        if tok.type == tokenize.COMMENT:
            continue
        out_tokens.append(tok)
    result = tokenize.untokenize(out_tokens)
    cleaned, prev_blank = ([], False)
    for line in result.splitlines():
        if line.strip() == '':
            if not prev_blank:
                cleaned.append('')
            prev_blank = True
        else:
            cleaned.append(line.rstrip())
            prev_blank = False
    return '\n'.join(cleaned) + '\n'

def strip_file(src: Path, dst: Path) -> None:
    source = src.read_text(encoding='utf-8')
    try:
        cleaned = remove_comments(remove_docstrings(source))
    except SyntaxError as e:
        print(f'  [SKIP] {src} — синтаксическая ошибка: {e}')
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(cleaned, encoding='utf-8')
    print(f'  {src} -> {dst}  ({len(source)} -> {len(cleaned)} символов)')

def process_dir(src_dir: Path, dst_dir: Path) -> None:
    files = sorted(src_dir.rglob('*.py'))
    if not files:
        print(f'В {src_dir} нет .py файлов.')
        return
    print(f'Найдено {len(files)} файлов в {src_dir}:')
    for f in files:
        rel = f.relative_to(src_dir)
        strip_file(f, dst_dir / rel)

def main() -> None:
    args = sys.argv[1:]
    if len(args) == 0:
        src = dst = Path.cwd()
    elif len(args) == 1:
        src = dst = Path(args[0]).resolve()
    elif len(args) == 2:
        src = Path(args[0]).resolve()
        dst = Path(args[1]).resolve()
    else:
        print(__doc__)
        sys.exit(1)
    if not src.is_dir():
        print(f'Не папка: {src}')
        sys.exit(1)
    process_dir(src, dst)
if __name__ == '__main__':
    main()
