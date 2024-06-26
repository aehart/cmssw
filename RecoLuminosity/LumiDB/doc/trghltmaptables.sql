/*==========================TRGHLTMAP======================*/
CREATE TABLE "TRGHLTMAP"
("HLTKEY" VARCHAR2(128) NOT NULL ENABLE,
"HLTPATHNAME" VARCHAR2(256) NOT NULL ENABLE,
"L1SEED" VARCHAR2(1024) NOT NULL ENABLE,
 CONSTRAINT "TRGHLTMAP_PK" PRIMARY KEY ("HLTKEY","HLTPATHNAME","L1SEED")
)
PARTITION BY HASH(HLTKEY)
PARTITIONS 4
STORE IN (CMS_LUMI_PROD_DATA,CMS_LUMI_PROD_DATA,CMS_LUMI_PROD_DATA,CMS_LUMI_PROD_DATA)
;
GRANT SELECT ON "TRGHLTMAP" TO PUBLIC;
GRANT SELECT,INSERT,DELETE,UPDATE ON "TRGHLTMAP" TO CMS_LUMI_WRITER;
