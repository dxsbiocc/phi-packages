import assert from 'node:assert/strict'
import { existsSync, readdirSync, readFileSync, statSync } from 'node:fs'
import { join, relative, resolve } from 'node:path'
import { fileURLToPath } from 'node:url'
import test from 'node:test'

const root = process.env.PHI_PACKAGES_TEST_ROOT ?? fileURLToPath(new URL('..', import.meta.url))
const skillRoot = join(root, 'resources/skills')
const scriptExtensions = new Set(['.py', '.R', '.r', '.sh', '.js', '.mjs', '.pl'])
const namePattern = /^[a-z0-9][a-z0-9-]{0,63}$/

function skillNames() {
  return readdirSync(skillRoot)
    .filter(name => statSync(join(skillRoot, name)).isDirectory())
    .sort()
}

function frontmatter(text, file) {
  assert.equal(text.startsWith('---\n'), true, `${file} must start with YAML frontmatter`)
  const end = text.indexOf('\n---\n', 4)
  assert.notEqual(end, -1, `${file} frontmatter is missing the closing ---`)
  return text.slice(4, end)
}

function rawField(fm, key) {
  const match = fm.match(new RegExp(`^${key}:\\s*(.*)$`, 'm'))
  return match ? match[1].trim() : undefined
}

function descriptionOf(fm) {
  const value = rawField(fm, 'description')
  assert.ok(value !== undefined, 'description is required')
  if (!['>', '>-', '|', '|-', '|+'].includes(value)) return value.replace(/^["']|["']$/g, '')
  const lines = []
  const start = fm.indexOf(`description: ${value}`)
  const rest = fm.slice(start).split('\n').slice(1)
  for (const line of rest) {
    if (line === '') {
      lines.push('')
      continue
    }
    if (!/^\s/.test(line)) break
    lines.push(line.trim())
  }
  return value.startsWith('|') ? lines.join('\n') : lines.join(' ').replace(/\s+/g, ' ').trim()
}

function hasScripts(dir) {
  const scripts = join(dir, 'scripts')
  if (!existsSync(scripts)) return false
  const stack = [scripts]
  while (stack.length > 0) {
    const current = stack.pop()
    for (const entry of readdirSync(current, { withFileTypes: true })) {
      const full = join(current, entry.name)
      if (entry.isDirectory()) stack.push(full)
      else if (scriptExtensions.has(entry.name.slice(entry.name.lastIndexOf('.')))) return true
    }
  }
  return false
}

function relativeLinks(text) {
  const links = []
  for (const match of text.matchAll(/\[[^\]]+\]\(([^)\s]+)\)/g)) {
    const target = match[1]
    if (target.startsWith('http://') || target.startsWith('https://') || target.startsWith('#')) continue
    links.push(target.split('#')[0])
  }
  return links.filter(Boolean)
}

test('README skill counts match standalone skill directories', () => {
  const count = skillNames().length
  assert.match(readFileSync(join(root, 'README.md'), 'utf8'), new RegExp(`\\| 独立 skills \\| ${count} \\|`))
  assert.match(readFileSync(join(root, 'README.en.md'), 'utf8'), new RegExp(`\\| Standalone skills \\| ${count} \\|`))
})

test('find-skills routing table names every standalone skill', () => {
  const routing = readFileSync(join(skillRoot, 'find-skills/references/routing.md'), 'utf8')
  for (const name of skillNames()) {
    assert.match(routing, new RegExp(`\`${name}\``), `routing.md is missing ${name}`)
  }
})

for (const name of skillNames()) {
  test(`${name} follows the skill contract`, () => {
    const dir = join(skillRoot, name)
    const file = join(dir, 'SKILL.md')
    const text = readFileSync(file, 'utf8')
    const fm = frontmatter(text, relative(root, file))
    const skillName = rawField(fm, 'name')
    assert.equal(skillName, name)
    assert.match(skillName, namePattern)
    assert.equal(skillName.includes('--'), false)
    assert.equal(skillName.endsWith('-'), false)

    const description = descriptionOf(fm)
    assert.equal(description.length >= 1 && description.length <= 1024, true, `description length ${description.length}`)

    const metadata = rawField(fm, 'metadata')
    if (metadata !== undefined) {
      assert.equal(metadata.startsWith('{'), true, 'metadata must be a one-line string map')
      const parsed = JSON.parse(metadata)
      assert.equal(typeof parsed, 'object')
      for (const value of Object.values(parsed)) assert.equal(typeof value, 'string')
    }

    const environmentDeclared = /^[ ]+environment:\s+\S+/m.test(fm)
    assert.equal(
      environmentDeclared,
      hasScripts(dir),
      environmentDeclared
        ? 'environment must not be declared when the skill has no scripts'
        : 'environment is required when the skill has scripts'
    )

    for (const target of relativeLinks(text)) {
      const resolved = resolve(dir, target)
      assert.equal(existsSync(resolved), true, `${name} links to missing ${target}`)
    }
  })
}
