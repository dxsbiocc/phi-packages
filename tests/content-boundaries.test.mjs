import assert from 'node:assert/strict'
import { existsSync, readFileSync } from 'node:fs'
import { join } from 'node:path'
import { fileURLToPath } from 'node:url'
import test from 'node:test'

const root = process.env.PHI_PACKAGES_TEST_ROOT ?? fileURLToPath(new URL('..', import.meta.url))

test('content source omits Phi-owned shared runtime and palette trees', () => {
  for (const relative of ['resources/palettes', 'resources/runtime']) {
    assert.equal(existsSync(join(root, relative)), false, `${relative} belongs to the Phi engine`)
  }
})

test('visualization retains the palette data its installed helpers read', () => {
  const skillRoot = join(root, 'resources/plugins/visualization/skills/omics-visualization')
  assert.ok(readFileSync(join(skillRoot, 'references/palettes.yaml'), 'utf8').length > 0)
  const catalog = JSON.parse(readFileSync(join(skillRoot, 'references/palettes/colors.json'), 'utf8'))
  assert.ok(catalog && typeof catalog === 'object')
})

const family = join(root, 'resources/wrappers/modules/local/differential-expression')
const dependencies = {
  deseq2: ['bioconductor-deseq2', 'bioconductor-limma', 'bioconductor-biocparallel', 'r-ashr'],
  edger: ['bioconductor-edger', 'bioconductor-limma'],
  limma: ['bioconductor-edger', 'bioconductor-limma'],
  qc: ['bioconductor-deseq2', 'r-ggplot2', 'r-pheatmap'],
  timeseries: ['bioconductor-variancepartition', 'bioconductor-edger', 'bioconductor-limma', 'bioconductor-biocparallel', 'r-lme4'],
  'visualization/pca': ['bioconductor-deseq2', 'r-ggplot2'],
  'visualization/heatmap': ['r-pheatmap'],
  'visualization/volcano': ['r-ggplot2']
}

test('retired differential-expression image source is absent', () => {
  assert.equal(existsSync(join(root, 'resources/wrappers/images/differential-expression-r')), false)
})

for (const [name, required] of Object.entries(dependencies)) {
  test(`${name} keeps executable code and declares its required R dependencies locally`, () => {
    const moduleRoot = join(family, name)
    const source = readFileSync(join(moduleRoot, 'main.nf'), 'utf8')
    assert.ok(!/phi\/differential-expression-r|images\/differential-expression-r/.test(source), 'retired image reference remains')
    const reference = source.match(/conda\s+["']\$\{moduleDir\}\/([^"']+)["']/)
    assert.ok(reference, 'module needs a package-relative Conda environment')
    assert.equal(reference[1].includes('..'), false, 'environment must stay with its module')
    const spec = readFileSync(join(moduleRoot, reference[1]), 'utf8')
    assert.match(spec, /^dependencies:\s*$/m)
    const packages = new Set([...spec.matchAll(/^\s*-\s*(?:[\w-]+::)?([\w-]+)(?:[=<>][^\n#]+)?\s*(?:#.*)?$/gm)].map(match => match[1]))
    for (const pkg of required) assert.ok(packages.has(pkg), `missing required dependency ${pkg}`)
    assert.ok(existsSync(join(moduleRoot, 'wrapper/wrapper.yaml')), 'keep the callable wrapper')
    assert.ok(existsSync(join(moduleRoot, 'wrapper/main.nf')), 'keep the wrapper entrypoint')
  })
}
