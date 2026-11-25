import assert from 'node:assert/strict'
import { existsSync, readdirSync, readFileSync, statSync } from 'node:fs'
import { dirname, join, resolve } from 'node:path'
import { fileURLToPath } from 'node:url'
import test from 'node:test'

const root = process.env.PHI_PACKAGES_TEST_ROOT ?? fileURLToPath(new URL('..', import.meta.url))

const ALLOWED_AGENT_KEYS = new Set([
  'name',
  'description',
  'tools',
  'model',
  'thinkingLevel',
  'skills',
  'environment',
  'visibility',
  'delegationMode',
  'delegation',
  'fallback',
])
const RESERVED_AGENT_KEYS = new Set(['spawns', 'outputSchema'])
const CONNECTOR_CATEGORIES = new Set(['生产力', '沟通协作', '设计创作', '健康与生命科学', '科研数据'])
const PLATFORMS = ['darwin-arm64', 'darwin-x64', 'linux-x64']
const RESERVED_PLUGIN_DIRS = ['mcp', 'wrappers', 'orchestrator']

function listDirs(dir) {
  return readdirSync(dir)
    .filter(n => n !== '.DS_Store' && statSync(join(dir, n)).isDirectory())
    .sort()
}

function extractFrontmatter(text, fileLabel) {
  assert.equal(text.startsWith('---\n'), true, `${fileLabel}: missing YAML frontmatter`)
  const end = text.indexOf('\n---\n', 4)
  assert.notEqual(end, -1, `${fileLabel}: unclosed YAML frontmatter`)
  return text.slice(4, end)
}

function topYamlKeys(yamlText) {
  return [...yamlText.matchAll(/^([A-Za-z_][A-Za-z0-9_-]*):/gm)].map(m => m[1])
}

function yamlScalar(yamlText, key) {
  const m = yamlText.match(new RegExp(`^${key}:\\s*(.*)$`, 'm'))
  if (!m) return undefined
  const raw = m[1].trim()
  if (!['>', '>-', '|', '|-', '|+'].includes(raw)) {
    return raw.replace(/^["']|["']$/g, '')
  }
  const idx = yamlText.indexOf(m[0])
  const rest = yamlText.slice(idx).split('\n').slice(1)
  const lines = []
  for (const line of rest) {
    if (line === '') {
      lines.push('')
      continue
    }
    if (!/^\s/.test(line)) break
    lines.push(line.trim())
  }
  return raw.startsWith('|') ? lines.join('\n') : lines.join(' ').replace(/\s+/g, ' ').trim()
}

test('all plugins, agents, and private environments satisfy docs/contracts/*', () => {
  const pluginsDir = join(root, 'resources/plugins')
  const plugins = listDirs(pluginsDir)
  const seenToolPrefixes = new Map()
  const seenAgentNames = new Map()
  const seenSkillNames = new Map()

  for (const s of listDirs(join(root, 'resources/skills'))) {
    seenSkillNames.set(s, `standalone:${s}`)
  }

  for (const plugin of plugins) {
    const pdir = join(pluginsDir, plugin)
    for (const reserved of RESERVED_PLUGIN_DIRS) {
      assert.equal(existsSync(join(pdir, reserved)), false, `${plugin}: reserved directory '${reserved}/' is forbidden`)
    }

    const pkgText = readFileSync(join(pdir, 'phi-package.yaml'), 'utf8')
    assert.equal(yamlScalar(pkgText, 'schemaVersion'), '1', `${plugin}: schemaVersion must be 1`)
    assert.equal(yamlScalar(pkgText, 'id'), plugin, `${plugin}: id must match directory`)
    assert.equal(yamlScalar(pkgText, 'type'), 'plugin', `${plugin}: type must be plugin`)
    assert.match(yamlScalar(pkgText, 'version') ?? '', /^\d+\.\d+\.\d+(-[0-9A-Za-z.-]+)?$/, `${plugin}: invalid semver`)

    const title = yamlScalar(pkgText, 'title') ?? ''
    assert.ok(title.length >= 1 && title.length <= 80, `${plugin}: title length ${title.length} outside 1..80`)
    const summary = yamlScalar(pkgText, 'summary') ?? ''
    assert.ok(summary.length >= 1 && summary.length <= 300, `${plugin}: summary length ${summary.length} outside 1..300`)

    const toolPrefix = yamlScalar(pkgText, 'toolPrefix') ?? ''
    assert.match(toolPrefix, /^[a-z][a-z0-9]{1,11}$/, `${plugin}: invalid toolPrefix '${toolPrefix}'`)
    assert.equal(seenToolPrefixes.has(toolPrefix), false, `duplicate toolPrefix '${toolPrefix}' between ${seenToolPrefixes.get(toolPrefix)} and ${plugin}`)
    seenToolPrefixes.set(toolPrefix, plugin)

    // Listed components vs filesystem
    const listedAgents = [...pkgText.matchAll(/^\s*-\s*(agents\/[A-Za-z0-9]+\.md)\s*$/gm)].map(m => m[1]).sort()
    const listedSkills = [...pkgText.matchAll(/^\s*-\s*(skills\/[a-z0-9-]+)\s*$/gm)].map(m => m[1]).sort()
    assert.ok(listedAgents.length + listedSkills.length >= 1, `${plugin}: must declare >=1 component`)

    const actualAgents = existsSync(join(pdir, 'agents'))
      ? readdirSync(join(pdir, 'agents'))
          .filter(f => f.endsWith('.md'))
          .map(f => `agents/${f}`)
          .sort()
      : []
    const actualSkills = existsSync(join(pdir, 'skills'))
      ? listDirs(join(pdir, 'skills')).map(s => `skills/${s}`).sort()
      : []
    assert.deepEqual(listedAgents, actualAgents, `${plugin}: components.agents mismatch with agents/ directory`)
    assert.deepEqual(listedSkills, actualSkills, `${plugin}: components.skills mismatch with skills/ directory`)

    // Environments & locks
    const declaredEnvs = new Set()
    const envDirRoot = join(pdir, 'environments')
    if (existsSync(envDirRoot)) {
      for (const envName of listDirs(envDirRoot)) {
        declaredEnvs.add(envName)
        assert.match(envName, /^[a-z][a-z0-9-]{0,62}$/, `${plugin}: invalid env name '${envName}'`)
        assert.match(pkgText, new RegExp(`^\\s{2}${envName}:\\s*$`, 'm'), `${plugin}: env '${envName}' missing in phi-package.yaml`)
        const specText = readFileSync(join(envDirRoot, envName, 'environment.yml'), 'utf8')
        assert.equal(yamlScalar(specText, 'name'), envName, `${plugin}/${envName}: environment.yml name mismatch`)
        for (const plat of PLATFORMS) {
          const lockPath = join(envDirRoot, envName, 'locks', `${plat}.txt`)
          assert.equal(existsSync(lockPath), true, `${plugin}/${envName}: missing lock ${plat}.txt`)
          const lines = readFileSync(lockPath, 'utf8')
            .split(/\r?\n/)
            .map(l => l.trim())
            .filter(l => l && !l.startsWith('#'))
          assert.equal(lines[0], '@EXPLICIT', `${plugin}/${envName}/${plat}: first non-comment line must be @EXPLICIT`)
          for (const urlLine of lines.slice(1)) {
            assert.match(urlLine, /^https:\/\/.+#[0-9a-f]{32}$/, `${plugin}/${envName}/${plat}: invalid lock URL '${urlLine}'`)
          }
        }
      }
    }

    // Validate every Agent file against docs/contracts/agent.md
    for (const relAgent of listedAgents) {
      const fileName = relAgent.slice('agents/'.length, -'.md'.length)
      const fm = extractFrontmatter(readFileSync(join(pdir, relAgent), 'utf8'), `${plugin}/${relAgent}`)
      const keys = topYamlKeys(fm)
      for (const k of keys) {
        assert.equal(RESERVED_AGENT_KEYS.has(k), false, `${plugin}/${relAgent}: reserved key '${k}' is rejected in v1`)
        assert.equal(ALLOWED_AGENT_KEYS.has(k), true, `${plugin}/${relAgent}: unknown or snake_case key '${k}'`)
      }
      const agentName = yamlScalar(fm, 'name')
      assert.equal(agentName, fileName, `${plugin}/${relAgent}: name '${agentName}' must equal filename '${fileName}'`)
      assert.match(agentName ?? '', /^[A-Z][A-Za-z0-9]*$/, `${plugin}/${relAgent}: invalid PascalCase name`)
      assert.equal(agentName.endsWith('Agent'), false, `${plugin}/${relAgent}: name must not end with 'Agent'`)
      assert.equal(seenAgentNames.has(agentName), false, `duplicate agent name '${agentName}'`)
      seenAgentNames.set(agentName, plugin)

      const desc = yamlScalar(fm, 'description') ?? ''
      assert.ok(desc.length >= 1 && desc.length <= 1024, `${plugin}/${relAgent}: description length ${desc.length} outside 1..1024`)
      assert.match(fm, /^tools:\s*(\n\s+-\s+\S+|\s*\[.+\])/m, `${plugin}/${relAgent}: non-empty tools list is required`)

      const vis = yamlScalar(fm, 'visibility')
      if (vis !== undefined) assert.equal(vis, 'entry', `${plugin}/${relAgent}: visibility '${vis}' is reserved`)

      const envRef = yamlScalar(fm, 'environment')
      if (envRef && envRef.startsWith('plugin:')) {
        const targetEnv = envRef.slice('plugin:'.length)
        assert.equal(declaredEnvs.has(targetEnv), true, `${plugin}/${relAgent}: undeclared plugin environment '${envRef}'`)
      }
    }

    // Validate every Skill script target, output schema existence, and private env reference
    const pluginAgentSet = new Set(listedAgents.map(a => a.slice('agents/'.length, -'.md'.length)))
    pluginAgentSet.add('main')
    for (const relSkill of listedSkills) {
      const sname = relSkill.slice('skills/'.length)
      assert.equal(seenSkillNames.has(sname), false, `duplicate skill name '${sname}'`)
      seenSkillNames.set(sname, plugin)

      const sdir = join(pdir, relSkill)
      const fm = extractFrontmatter(readFileSync(join(sdir, 'SKILL.md'), 'utf8'), `${plugin}/${relSkill}`)
      const sEnv = (fm.match(/^\s+environment:\s+(\S+)/m) ?? [])[1]
      if (sEnv && sEnv.startsWith('plugin:')) {
        assert.equal(declaredEnvs.has(sEnv.slice('plugin:'.length)), true, `${plugin}/${relSkill}: undeclared environment '${sEnv}'`)
      }
      for (const m of fm.matchAll(/^\s*run:\s*\[([^\]]+)\]/gm)) {
        const parts = m[1].split(',').map(x => x.trim().replace(/^["']|["']$/g, ''))
        for (const arg of parts.slice(1)) {
          if (arg.startsWith('./')) {
            assert.equal(existsSync(join(sdir, arg)), true, `${plugin}/${relSkill}: missing script '${arg}'`)
          }
        }
      }
      const toolCount = (fm.match(/^    - name:\s+\S+/gm) ?? []).length
      const outMatches = [...fm.matchAll(/^      output:\s+(\S+)/gm)]
      assert.equal(outMatches.length, toolCount, `${plugin}/${relSkill}: ${outMatches.length}/${toolCount} script tools declare an output JSON schema`)
      for (const m of outMatches) {
        assert.equal(existsSync(join(sdir, m[1])), true, `${plugin}/${relSkill}: missing output schema '${m[1]}'`)
        JSON.parse(readFileSync(join(sdir, m[1]), 'utf8'))
      }
    }
  }
})

test('all MCP connectors satisfy docs/contracts/connector.md', () => {
  const connRoot = join(root, 'resources/connectors')
  for (const cid of listDirs(connRoot)) {
    const cdir = join(connRoot, cid)
    const txt = readFileSync(join(cdir, 'phi-package.yaml'), 'utf8')
    assert.equal(yamlScalar(txt, 'schemaVersion'), '1', `${cid}: schemaVersion must be 1`)
    assert.equal(yamlScalar(txt, 'id'), cid, `${cid}: id must equal directory`)
    assert.equal(yamlScalar(txt, 'type'), 'mcp', `${cid}: type must be mcp`)
    const title = yamlScalar(txt, 'title') ?? ''
    const summary = yamlScalar(txt, 'summary') ?? ''
    assert.ok(title.length >= 1 && title.length <= 80, `${cid}: invalid title length`)
    assert.ok(summary.length >= 1 && summary.length <= 300, `${cid}: invalid summary length`)

    const transport = (txt.match(/^\s+transport:\s*(\S+)/m) ?? [])[1]
    assert.ok(transport === 'http' || transport === 'stdio', `${cid}: invalid transport '${transport}'`)
    const cat = (txt.match(/^\s+category:\s*(.+)$/m) ?? [])[1]?.trim().replace(/^["']|["']$/g, '')
    assert.equal(CONNECTOR_CATEGORIES.has(cat), true, `${cid}: invalid category '${cat}'`)
    assert.equal(/^\s+secrets:/m.test(txt), false, `${cid}: connector.secrets is reserved in v1`)
    if (transport === 'http') {
      const url = (txt.match(/^\s+url:\s*(\S+)/m) ?? [])[1]?.replace(/^["']|["']$/g, '')
      assert.match(url ?? '', /^https:\/\//, `${cid}: http connector requires https:// url`)
      const auth = (txt.match(/^\s+auth:\s*(\S+)/m) ?? [])[1]
      assert.ok(['none', 'oauth', 'header'].includes(auth), `${cid}: invalid auth '${auth}'`)
    }
  }
})

test('all 2202 wrapper adapters keep required files and all Nextflow relative includes resolve', () => {
  const wrapRoot = join(root, 'resources/wrappers')
  let adapterCount = 0
  let includeCount = 0
  const brokenIncludes = []

  const stack = [wrapRoot]
  while (stack.length > 0) {
    const cur = stack.pop()
    for (const e of readdirSync(cur, { withFileTypes: true })) {
      if (e.name === '.DS_Store') continue
      const full = join(cur, e.name)
      if (e.isDirectory()) {
        stack.push(full)
      } else if (e.name === 'wrapper.yaml' && cur.endsWith('/wrapper')) {
        adapterCount += 1
        assert.equal(existsSync(join(cur, 'main.nf')), true, `${cur}: missing main.nf`)
        assert.equal(existsSync(join(cur, 'params.json')), true, `${cur}: missing params.json`)
        JSON.parse(readFileSync(join(cur, 'params.json'), 'utf8'))
      } else if (e.name.endsWith('.nf')) {
        const nfText = readFileSync(full, 'utf8')
        for (const m of nfText.matchAll(/^\s*include\s*\{[^}]+\}\s*from\s*["']([^"']+)["']/gm)) {
          const target = m[1]
          if (target.startsWith('plugin/')) continue
          includeCount += 1
          const resolved = resolve(dirname(full), target)
          if (!existsSync(resolved) && !existsSync(`${resolved}.nf`) && !existsSync(join(resolved, 'main.nf'))) {
            brokenIncludes.push(`${full} -> ${target}`)
          }
        }
      }
    }
  }

  assert.ok(adapterCount >= 2200, `unexpected wrapper count ${adapterCount}`)
  assert.ok(includeCount >= 4000, `unexpected include count ${includeCount}`)
  assert.deepEqual(brokenIncludes, [], `broken Nextflow includes:\n${brokenIncludes.slice(0, 10).join('\n')}`)
})
