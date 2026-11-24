#!/usr/bin/env python3
"""Standalone YAML adapter to the shared gsim R/native workflow; no R runtime."""
import argparse,json,math,os,re,subprocess,sys,tempfile,time
from pathlib import Path
import yaml
class Loader(yaml.SafeLoader):pass
def mapping(loader,node):
 out={}
 for k,v in node.value:
  key=loader.construct_object(k)
  if not isinstance(key,str) or key in out:raise ValueError('YAML keys must be unique strings')
  out[key]=loader.construct_object(v,deep=True)
 return out
Loader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,mapping)
def wire(path):
 path=Path(path).resolve(strict=True);text=path.read_text(encoding='utf-8-sig')
 if any(isinstance(t,(yaml.tokens.AnchorToken,yaml.tokens.AliasToken,yaml.tokens.TagToken)) for t in yaml.scan(text)):raise ValueError('YAML anchors, aliases and tags are unsupported')
 doc=yaml.load(text,Loader=Loader)
 keys={'schema','inputs','traits','seed','genotyped','phenotyped','records_per_individual','visits','missing_rate','residual_covariance','maximum_memory_bytes','export_truth_phases'}
 if isinstance(doc,dict) and doc.get('schema')=='gsim.gaussian_model/v1':
  if set(doc)-{'schema','inputs','seed','maximum_memory_bytes'} or not isinstance(doc.get('inputs'),dict) or set(doc['inputs'])!={'specification'}:raise ValueError('prepared Gaussian recipe requires inputs.specification and optional seed/memory controls')
  value=doc['inputs']['specification']
  if not isinstance(value,str) or not value:raise ValueError('inputs.specification must be a nonempty path')
  spec=(path.parent/value).resolve(strict=True)
  if not spec.is_file():raise ValueError('inputs.specification must be a regular file')
  out='GSIM_GAUSSIAN 1\nspecification '+json.dumps(spec.as_posix(),ensure_ascii=False)+'\n'
  for key in ('seed','maximum_memory_bytes'):
   if key in doc:
    value=doc[key]
    if isinstance(value,bool) or not isinstance(value,int) or not 0<=value<=2**53-3:raise ValueError('invalid Gaussian integer control: '+key)
    out+=key+' '+json.dumps(str(value))+'\n'
  return out
 if not isinstance(doc,dict) or set(doc)-keys or doc.get('schema')!='gsim.simulation/v1':raise ValueError('invalid simulation schema or unknown keys')
 if 'records_per_individual' in doc and 'visits' in doc:raise ValueError('specify records_per_individual or legacy visits, not both')
 inputs=doc.get('inputs');traits=doc.get('traits')
 if not isinstance(inputs,dict) or set(inputs)!={'pedigree','variants','founders','effects'}:raise ValueError('inputs must name pedigree, variants, founders and effects')
 if not isinstance(traits,list) or not traits or any(not isinstance(t,str) for t in traits):raise ValueError('traits must be a nonempty string list')
 if any(not re.fullmatch(r'[A-Za-z0-9_.-]+',t) or t=='0' for t in traits):raise ValueError('invalid trait identity')
 c={}
 for k,v in inputs.items():
  if not isinstance(v,str) or not v:raise ValueError('input paths must be nonempty strings')
  p=(path.parent/v).resolve(strict=True)
  if not p.is_file():raise ValueError('input must be a regular file')
  c[k]=p.as_posix()
 if 'export_truth_phases' in doc:
  if not isinstance(doc['export_truth_phases'],bool):raise ValueError('export_truth_phases must be boolean')
  c['export_truth_phases']='1' if doc['export_truth_phases'] else '0'
 c['traits']=' '.join(traits)
 for k in ('seed','genotyped','phenotyped','records_per_individual','visits','maximum_memory_bytes'):
  if k in doc:
   v=doc[k]
   if isinstance(v,bool) or not isinstance(v,int):raise ValueError(k+' must be an integer')
   c[k]=str(v)
 if 'missing_rate' in doc:
  v=doc['missing_rate']
  if isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v):raise ValueError('missing_rate must be finite numeric')
  c['missing_rate']=format(v,'.17g')
 if 'residual_covariance' in doc:
  a=doc['residual_covariance'];q=len(traits)
  if not isinstance(a,list) or len(a)!=q or any(not isinstance(row,list) or len(row)!=q for row in a):raise ValueError('residual covariance must be a trait-ordered square matrix')
  if any(isinstance(x,bool) or not isinstance(x,(int,float)) or not math.isfinite(x) for row in a for x in row):raise ValueError('residual covariance must be finite numeric')
  c['residual_covariance']=' '.join(format(x,'.17g') for row in a for x in row)
 return 'GSIM_SIMULATION 1\n'+''.join(k+' '+json.dumps(v,ensure_ascii=False)+'\n' for k,v in c.items())
def main():
 p=argparse.ArgumentParser(prog='gsim');cmd=p.add_subparsers(dest='command',required=True);x=cmd.add_parser('simulate');x.add_argument('recipe',type=Path);x.add_argument('--output',required=True,type=Path);x.add_argument('--exe',type=Path);x.add_argument('--model-driver',type=Path);x.add_argument('--model-preparer',type=Path)
 a=p.parse_args();started=time.perf_counter()
 document=yaml.load(a.recipe.read_text(encoding='utf-8-sig'),Loader=Loader)
 full_gaussian=isinstance(document,dict) and document.get('schema')=='gsim.gaussian_model/v1' and 'model' in document
 recipe=None if full_gaussian else wire(a.recipe)
 exe=a.exe
 if exe is None:
  for suffix in ('','.exe'):
   candidate=Path(__file__).resolve().parents[2]/'bin'/('gsim_native'+suffix)
   if candidate.is_file():exe=candidate;break
 if exe is None:raise ValueError('gsim_native not installed; supply --exe or install the standalone build')
 target=a.output.resolve()
 if target.exists():raise ValueError('preserve existing simulation output')
 if full_gaussian:
  allowed={'schema','model','design_inputs','seed','maximum_memory_bytes','maximum_dimension'}
  text=a.recipe.read_text(encoding='utf-8-sig')
  if any(isinstance(t,(yaml.tokens.AnchorToken,yaml.tokens.AliasToken,yaml.tokens.TagToken)) for t in yaml.scan(text)):raise ValueError('YAML anchors, aliases and tags are unsupported')
  if set(document)-allowed or not isinstance(document.get('model'),str) or not isinstance(document.get('design_inputs'),dict) or not document['design_inputs']:raise ValueError('invalid Gaussian model simulation recipe')
  for key in ('seed','maximum_memory_bytes','maximum_dimension'):
   if key in document and (isinstance(document[key],bool) or not isinstance(document[key],int) or not 0<=document[key]<=2**53-3):raise ValueError('invalid Gaussian integer control: '+key)
  root=Path(__file__).resolve().parents[2];preparer=a.model_preparer or root/'share/gsuite/model_simulation.py'
  driver=a.model_driver
  if driver is None:
   for suffix in ('','.exe'):
    candidate=root/'bin'/('gsuite_model_driver'+suffix)
    if candidate.is_file():driver=candidate;break
  if not preparer.is_file() or driver is None:raise ValueError('Gaussian model generation requires matching gsuite preparation and native driver; install the matched bundle or supply --model-preparer and --model-driver')
  base=a.recipe.resolve(strict=True).parent;model=(base/document['model']).resolve(strict=True)
  design=[]
  for key,value in document['design_inputs'].items():
   if not isinstance(key,str) or not isinstance(value,str) or not value:raise ValueError('invalid Gaussian design input')
   design+=['--design-input',key+'='+str((base/value).resolve(strict=True))]
  prepared=target.with_name(target.name+'.preparation')
  args=[sys.executable,str(preparer),'prepare','--model',str(model),*design,'--output',str(prepared),'--exe',str(driver)]
  if 'maximum_dimension' in document:args+=['--maximum-dimension',str(document['maximum_dimension'])]
  subprocess.run(args,check=True)
  recipe='GSIM_GAUSSIAN 1\nspecification '+json.dumps((prepared/'specification.txt').as_posix(),ensure_ascii=False)+'\n'
  for key in ('seed','maximum_memory_bytes'):
   if key in document:recipe+=key+' '+json.dumps(str(document[key]))+'\n'
 with tempfile.TemporaryDirectory(prefix='gsim-recipe-') as temp:
  config=Path(temp)/'recipe.txt';config.write_text(recipe,encoding='utf8',newline='\n')
  subprocess.run([str(exe.resolve(strict=True)),'--wire',str(config),str(target)],check=True)
 if full_gaussian:subprocess.run([sys.executable,str(preparer),'collect','--model',str(model),*design,'--prepared',str(prepared),'--output',str(target)],check=True)
 print(json.dumps(dict(accepted=True,output=str(target),simulation_seconds=time.perf_counter()-started)))
if __name__=='__main__':
 try:main()
 except (ValueError,OSError,subprocess.CalledProcessError,yaml.YAMLError) as e:raise SystemExit('gsim: '+str(e))
