for dataset in holatr_CI_SkillingsH dasalb_CR_Bernardi_Ramon abuabd_CB_coleman_newH abuvai_CB_coleman_newH acanigF_CB_eble acanigR_CB_DiBattista_new acaoli_CB_gaither acapla_CR_timmers anon3_CI_blank anon4_CI_blank anon8_CI_blank anon9_CR_blank carmel_ATP_Santos celexa_CI_Bird celsan_CI_Bird celtal_CI_Bird chafre_CB_craig_newH chalun_CB_Szabo chamil_CB_craigH chamul_CB_craigH ctestr_CB_eble epique_CR_Andrews etemar_CB_andrews hetmam_CI_bollick holwhi_CI_SkillingsH lutkas_CB_gaither mulfla_CB mulvan_CB myrber_CB_CraigH opheri_16Snew ophpic_16S_Skillings panmar_CI_Iacchei panpen_CI_Iacchei parmul_CB_szabo prifil_CB_gaither squmit_CR_dalyengel stefas_CR_Ramon stelon_CR_Andrews triobe_CR_Whitney zebfla_CB_Eble_nolump
  do
  echo $dataset
  cd /home/edc42/old_lus/Migrate_runs_round6/$dataset
  ../../submitter_UH.sh
  cd /home/edc42/old_lus/Migrate_runs_round6
done